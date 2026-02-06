Давай посмотрю твою текущую реализацию метода Ньютона в notebook'е, чтобы понять что именно ты используешь и предложить классическую схему.Файл очень большой. Давай найду конкретную реализацию метода Ньютона:Теперь поищу C-реализацию в notebook'е:Поищу существующую реализацию newton_step в файле mendrive_newton.c:Теперь всё понятно. Твой "4-step Newton" — это разделение комплексной функции на 4 компоненты:
1. Re(f/df_dkz) → коррекция Re(z) (kz)
2. Im(f/df_dkz) → коррекция Re(z) (kz)
3. Re(f/df_dsz) → коррекция Im(z) (sz)
4. Im(f/df_dsz) → коррекция Im(z) (sz)

Это своеобразный "coordinate descent" с адаптивным размером шага, а не стандартный комплексный Ньютон.

Классический комплексный метод Ньютона решает систему 2×2:
```
[∂f_re/∂x  ∂f_re/∂y] [Δx]   [f_re]
[∂f_im/∂x  ∂f_im/∂y] [Δy] = [f_im]
```

Создам классическую реализацию с поддержкой разных типов точности:Теперь создам адаптер для интеграции с твоей библиотекой mendrive_det:Теперь создам пример использования и тестовую программу:Скомпилируем и протестируем:Отлично работает! Теперь тест с double:Всё работает. Скопирую файлы в outputs:Вот три файла:

**Сравнение твоего метода и классического:**

| Аспект | Твой метод (4-step) | Классический Newton |
|--------|---------------------|---------------------|
| Шаг | 4 раздельных шага по координатам | 1 совместный шаг |
| Якобиан | Использует f/f' напрямую | Решает систему 2×2 |
| Линейная система | Не решает | J·Δz = -f(z) |
| Сходимость | Coordinate descent | Квадратичная вблизи корня |
| Устойчивость | Адаптивный размер по каждой координате | Damping/backtracking |

**Использование:**

```c
// Компиляция с твоей библиотекой:
gcc -shared -fPIC -O3 -DUSE_LONG_DOUBLE -DKY \
    mendrive_det.c \
    newton_mendrive_classic.c \
    ... \
    -o libmendrive.so -lm

// В Python:
lib.newton_classic_solve.argtypes = [
    POINTER(c_longdouble), POINTER(c_longdouble),
    POINTER(c_longdouble), POINTER(c_int)
]
lib.newton_classic_solve.restype = c_int

kz = c_longdouble(kz0)
sz = c_longdouble(sz0)
f_abs = c_longdouble()
niter = c_int()
status = lib.newton_classic_solve(byref(kz), byref(sz), byref(f_abs), byref(niter))
```

Флаги точности: `-DUSE_LONG_DOUBLE` или `-DUSE_FLOAT128` (с `-lquadmath`).