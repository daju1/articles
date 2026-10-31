#!/usr/bin/env python3
"""
nb_run.py — headless-выполнение Jupyter-ноутбука (SageMath и любые другие ядра).

* сохраняет результат после каждой ячейки, атомарно (tmp-файл + os.replace);
* пишет вывод ячеек в лог в реальном времени, с метками времени и длительностью;
* при ошибке в ячейке, смерти ядра, SIGTERM или SIGINT сохраняет частичный
  результат и корректно останавливает ядро (1-й сигнал — interrupt ядра,
  2-й — жёсткая остановка);
* защищён от повторного запуска блокировкой fcntl.flock; в lock-файле лежит PID;
* не пишет в ячейки служебные метки времени (metadata.execution), а уже
  имеющиеся в исходном ноутбуке удаляет; вернуть их можно ключом --keep-timing;
* перед расчётом очищает выводы и номера выполнения всех ячеек с кодом, поэтому
  в результате (X_running.ipynb) есть только то, что посчитано в этом запуске;
  невыполненные ячейки остаются пустыми, а исходный ноутбук не меняется
  (если нет --replace).

Коды выхода:
  0   успех
  1   ошибка в ячейке
  2   ядро упало или не запустилось
  3   прочая ошибка
  4   не удалось сохранить результат
  75  расчёт для этого ноутбука уже идёт
  128+N  остановлен сигналом N (143 = SIGTERM, 130 = SIGINT)
"""
import argparse
import asyncio
import datetime as dt
import fcntl
import os
import re
import signal
import sys
import tempfile
import time
import traceback

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError, DeadKernelError
from nbclient.util import ensure_async

ANSI_RE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
_log_state = {"at_line_start": True}


# ----------------------------------------------------------------- логирование
def log(msg: str) -> None:
    if not _log_state["at_line_start"]:
        print(flush=True)
    stamp = dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{stamp}] {msg}", flush=True)
    _log_state["at_line_start"] = True


def raw(text: str) -> None:
    if not text:
        return
    sys.stdout.write(text)
    sys.stdout.flush()
    _log_state["at_line_start"] = text.endswith("\n")


def fmt_dur(seconds: float) -> str:
    s = int(round(seconds))
    h, rem = divmod(s, 3600)
    m, s = divmod(rem, 60)
    return f"{h}:{m:02d}:{s:02d}" if h else f"{m}:{s:02d}"


class StreamingClient(NotebookClient):
    """Печатает вывод ячейки в лог по мере поступления, а не после её завершения."""

    def output(self, outs, msg, display_id, cell_index):
        msg_type = msg.get("msg_type")
        content = msg.get("content", {})
        if msg_type == "stream":
            raw(content.get("text", ""))
        elif msg_type in ("execute_result", "display_data"):
            data = content.get("data", {})
            if "text/plain" in data:
                raw(str(data["text/plain"]) + "\n")
            elif any(k.startswith("image/") for k in data):
                raw("[изображение]\n")
        elif msg_type == "error":
            log(f"[ОШИБКА В ЯЧЕЙКЕ {cell_index + 1}] "
                f"{content.get('ename')}: {content.get('evalue')}")
            for line in content.get("traceback", []):
                raw(ANSI_RE.sub("", line) + "\n")
        return super().output(outs, msg, display_id, cell_index)


# ----------------------------------------------------------------- файлы
def strip_execution(nb) -> int:
    """Удаляет metadata.execution (метки времени nbclient) из всех ячеек.

    Возвращает число ячеек, из которых поле было удалено.
    """
    removed = 0
    for cell in nb.cells:
        if "execution" in cell.get("metadata", {}):
            del cell["metadata"]["execution"]
            removed += 1
    return removed


def strip_outputs(nb) -> int:
    """Очищает выводы и execution_count у всех ячеек с кодом.

    Результат запуска не должен содержать выводов, оставшихся от прежних
    расчётов. Возвращает число ячеек, у которых были выводы или номер выполнения.
    """
    cleared = 0
    for cell in nb.cells:
        if cell.cell_type != "code":
            continue
        if cell.get("outputs") or cell.get("execution_count") is not None:
            cleared += 1
        cell["outputs"] = []
        cell["execution_count"] = None
    return cleared


def atomic_write(nb, path: str, mode: int) -> None:
    """Пишет ноутбук во временный файл рядом и атомарно подменяет им path."""
    directory = os.path.dirname(path)
    fd, tmp = tempfile.mkstemp(dir=directory,
                               prefix="." + os.path.basename(path) + ".",
                               suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as f:
            nbformat.write(nb, f)
        os.chmod(tmp, mode)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def acquire_lock(path: str):
    """Эксклюзивная блокировка; файл держим открытым до конца процесса."""
    f = open(path, "a+")
    try:
        fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except OSError:
        f.seek(0)
        other = f.read().strip() or "?"
        log(f"Расчёт уже идёт (PID {other}), блокировка: {path}")
        sys.exit(75)
    f.seek(0)
    f.truncate()
    f.write(str(os.getpid()))
    f.flush()
    return f


# ----------------------------------------------------------------- выполнение
async def run(client: StreamingClient, nb, args, state: dict) -> None:
    loop = asyncio.get_running_loop()
    task = asyncio.current_task()

    # Сигналы (SIGHUP намеренно не трогаем — его игнорирование от nohup сохраняется):
    #   во время запуска ядра — сразу остановка;
    #   1-й во время расчёта — interrupt ядра (как кнопка Interrupt): текущая ячейка
    #       получает KeyboardInterrupt, её вывод до этого момента сохраняется;
    #   2-й — жёсткая остановка (отмена задачи, ядро гасится).
    def on_signal(sig):
        name = signal.Signals(sig).name
        if state["signal"] is None:
            state["signal"] = sig
            if state["ready"]:
                log(f"Получен {name}: прерываю ядро (повторный сигнал — жёсткая остановка)")
                asyncio.ensure_future(ensure_async(client.km.interrupt_kernel()))
                return
            log(f"Получен {name} во время запуска ядра — остановка")
        else:
            log(f"Повторный {name}: жёсткая остановка")
        task.cancel()

    def install_handlers():
        for sig in (signal.SIGTERM, signal.SIGINT):
            loop.add_signal_handler(sig, on_signal, sig)

    install_handlers()

    # Ядро запускаем сами, под своими обработчиками: async_setup_kernel ставит
    # собственные (только гасят ядро, расчёт не останавливают), и на время долгого
    # старта Sage сигнал ушёл бы им.
    client.km = client.create_kernel_manager()
    try:
        await client.async_start_new_kernel()
        await client.async_start_new_kernel_client()
    except BaseException:
        await client._async_cleanup_kernel()
        raise

    async with client.async_setup_kernel():
        install_handlers()   # вернуть свои поверх только что поставленных nbclient
        state["ready"] = True

        log(f"Ядро '{client.kernel_name}' запущено.")
        try:
            reply = await client.async_wait_for_reply(client.kc.kernel_info())
            nb.metadata["language_info"] = reply["content"]["language_info"]
        except Exception:
            pass

        code_cells = [(i, c) for i, c in enumerate(nb.cells) if c.cell_type == "code"]
        total = len(code_cells)
        last_save = time.monotonic()

        for n, (idx, cell) in enumerate(code_cells, 1):
            if state["signal"] is not None:
                break
            state["cell"] = idx
            log(f">>> Ячейка {idx + 1} ({n}/{total})")
            if not args.no_source and cell.source.strip():
                for line in cell.source.splitlines():
                    raw(f"  | {line}\n")
                raw("  ---\n")

            t0 = time.monotonic()
            await client.async_execute_cell(cell, idx)
            log(f"<<< Ячейка {idx + 1} выполнена за {fmt_dur(time.monotonic() - t0)}")

            if time.monotonic() - last_save >= args.save_interval:
                atomic_write(nb, args.output, state["mode"])
                last_save = time.monotonic()


def parse_args():
    p = argparse.ArgumentParser(description="Headless-выполнение Jupyter-ноутбука")
    p.add_argument("input", help="исходный ноутбук (не изменяется, если нет --replace)")
    p.add_argument("output", help="куда писать результат (например, X_running.ipynb)")
    p.add_argument("--kernel", help="имя ядра; по умолчанию — из метаданных ноутбука")
    p.add_argument("--startup-timeout", type=int, default=180,
                   help="секунд на запуск ядра (Sage стартует долго), по умолчанию 180")
    p.add_argument("--cell-timeout", type=int, default=-1,
                   help="таймаут ячейки в секундах; -1 = без ограничения (по умолчанию)")
    p.add_argument("--save-interval", type=float, default=0.0,
                   help="не сохранять чаще, чем раз в N секунд (0 = после каждой ячейки)")
    p.add_argument("--no-source", action="store_true",
                   help="не печатать исходный код ячеек в лог")
    p.add_argument("--keep-timing", action="store_true",
                   help="записывать в ячейки метки времени выполнения (metadata.execution); "
                        "по умолчанию они не пишутся, а старые удаляются")
    p.add_argument("--lock", help="lock-файл; по умолчанию .<output>.lock рядом с output")
    p.add_argument("--replace", action="store_true",
                   help="после успешного завершения заменить input результатом")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.input = os.path.abspath(args.input)
    args.output = os.path.abspath(args.output)
    lock_path = args.lock or os.path.join(
        os.path.dirname(args.output), "." + os.path.basename(args.output) + ".lock")
    _lock = acquire_lock(lock_path)  # noqa: F841 — держим ссылку до конца процесса

    with open(args.input, encoding="utf-8") as f:
        nb = nbformat.read(f, as_version=4)

    # Метки времени, оставшиеся от прежних запусков, убираем до расчёта: иначе у
    # невыполненных ячеек они сохранились бы в результате.
    removed = 0 if args.keep_timing else strip_execution(nb)

    # Выводы прежних запусков тоже убираем: в результате должно остаться только
    # то, что прошло расчёт в этом запуске.
    cleared = strip_outputs(nb)

    kernel = args.kernel or nb.metadata.get("kernelspec", {}).get("name") or "python3"
    state = {"cell": None, "signal": None, "ready": False,
             "mode": os.stat(args.input).st_mode & 0o777}

    client = StreamingClient(
        nb,
        kernel_name=kernel,
        timeout=None if args.cell_timeout < 0 else args.cell_timeout,
        startup_timeout=args.startup_timeout,
        record_timing=args.keep_timing,   # False: nbclient не пишет metadata.execution
        resources={"metadata": {"path": os.path.dirname(args.input)}},  # cwd ядра
        extra_arguments=["--IPKernelApp.capture_fd_output=False"],  # вывод C-библиотек (Giac) не через трубу
    )

    log(f"Старт: {args.input}")
    log(f"Результат: {args.output}; ядро: {kernel}; PID {os.getpid()}")
    log(f"Очищены старые выводы в ячейках с кодом: {cleared}")
    if removed:
        log(f"Удалены старые метки времени (metadata.execution) в ячейках: {removed}")
    t_start = time.monotonic()
    code, status = 0, "успешно"

    try:
        asyncio.run(run(client, nb, args, state))
    except CellExecutionError:
        code, status = 1, "ошибка в ячейке"          # traceback уже выведен в лог
    except DeadKernelError as e:
        code, status = 2, f"ядро умерло ({e})"
    except (asyncio.CancelledError, KeyboardInterrupt):
        sig = state["signal"] or signal.SIGINT
        code, status = 128 + int(sig), f"остановлен сигналом {signal.Signals(sig).name}"
    except Exception as e:
        code, status = 3, f"{type(e).__name__}: {e}"
        traceback.print_exc()

    # Остановка по сигналу может проявиться как CellExecutionError (KeyboardInterrupt
    # в ячейке) или DeadKernelError (nbclient так оформляет отмену) — итог один.
    if state["signal"] is not None:
        sig = state["signal"]
        code, status = 128 + int(sig), f"остановлен сигналом {signal.Signals(sig).name}"

    if not args.keep_timing:
        strip_execution(nb)   # страховка: в результате меток времени быть не должно

    try:
        atomic_write(nb, args.output, state["mode"])
    except Exception as e:
        log(f"Не удалось сохранить результат: {e}")
        code = code or 4

    where = "" if state["cell"] is None or code == 0 else f", последняя ячейка {state['cell'] + 1}"
    log(f"Итог: {status}{where}; время {fmt_dur(time.monotonic() - t_start)}; код {code}")

    if code == 0 and args.replace:
        os.replace(args.output, args.input)
        log(f"Исходный ноутбук заменён результатом: {args.input}")

    sys.exit(code)


if __name__ == "__main__":
    main()
