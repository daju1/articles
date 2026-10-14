import nbformat
from nbclient import NotebookClient
import sys
import asyncio

input_path = sys.argv[1]
output_path = sys.argv[2]

# Читаем исходный ноутбук
with open(input_path, encoding='utf-8') as f:
    nb = nbformat.read(f, as_version=4)

# Получаем имя ядра из метаданных ноутбука
kernel_name = nb.metadata.get('kernelspec', {}).get('name', 'python3')
print(f"Используется ядро: {kernel_name}", flush=True)

# Инициализируем клиент
client = NotebookClient(
    nb, 
    timeout=-1,
    kernel_name=kernel_name,
    startup_timeout=180  # SageMath тяжелое, даём 3 минуты на запуск
)

last_executed_index = -1

def log_cell_input(cell, cell_index):
    """Выводит исходный код ячейки в лог перед выполнением."""
    source = cell.get('source', '')

    if not source.strip():
        print(f"  [Пустая ячейка]", flush=True)
        return

    print(f"  --- Исходный код ячейки {cell_index + 1} ---", flush=True)

    # Выводим код с отступами для читаемости
    for line in source.split('\n'):
        print(f"  | {line}", flush=True)

    print(f"  --- Конец кода ---", flush=True)

def log_cell_output(cell):
    """Выводит текстовые результаты ячейки, игнорируя изображения."""
    for output in cell.get('outputs', []):
        out_type = output.get('output_type')

        if out_type == 'stream':
            # Обычный вывод (print) или ошибки stderr
            # end='' используется, так как в text уже есть переносы строк
            print(output.get('text', ''), end='', flush=True)
 
        elif out_type == 'error':
            # Вывод traceback при ошибке
            print(f"\n[ОШИБКА В ЯЧЕЙКЕ] {output.get('ename')}: {output.get('evalue')}", flush=True)
            for line in output.get('traceback', []):
                print(line, flush=True)
   
        elif out_type in ('execute_result', 'display_data'):
            data = output.get('data', {})
            # Если есть текстовое представление, выводим его.
            # Это автоматически игнорирует ячейки, которые вернули ТОЛЬКО картинку (image/png)
            if 'text/plain' in data:
                print(data['text/plain'], flush=True)

async def run_notebook():
    global last_executed_index

    async with client.async_setup_kernel():
        print("Ядро успешно запущено.", flush=True)

        for index, cell in enumerate(nb.cells):
            if cell.cell_type == 'code':
                last_executed_index = index
                print(f"\n>>> Выполняется ячейка {index + 1}...", flush=True)

                # 1. Логируем исходный код ячейки
                log_cell_input(cell, index)

                # Выполняем ячейку
                await client.async_execute_cell(cell, index)
                
                # 1. Логируем текстовый вывод ячейки
                log_cell_output(cell)
                
                # 2. Сохраняем ноутбук после каждой выполненной ячейки
                with open(output_path, 'w', encoding='utf-8') as f:
                    nbformat.write(nb, f)
                print(f"<<< Ячейка {index + 1} выполнена и сохранена.\n", flush=True)
                
        print("Ноутбук успешно выполнен до конца.", flush=True)

try:
    asyncio.run(run_notebook())

except Exception as e:
    print(f"\nКРИТИЧЕСКАЯ ОШИБКА при выполнении ячейки {last_executed_index + 1}: {e}", flush=True)
    print("Сохраняем промежуточный результат перед выходом...", flush=True)
    
    with open(output_path, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)
        
    sys.exit(1)

finally:
    # ⚠️ Обязательно останавливаем ядро, иначе останутся висящие процессы sage
    try:
        client.cleanup_kernel()
        print("Ядро успешно остановлено.", flush=True)
        # Очистка уже происходит при выходе из async with, но явный вызов безопасен
        pass
    except Exception:
        pass