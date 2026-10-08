from ..variables.common import x, y, z, b
from ..common import get_use_phase_y

def avg_over_y(expr, y_min=-b, y_max=b):
    if get_use_phase_y():
        """Усреднение выражения по y с нормировкой на длину интервала."""
        # from sage.calculus.functional import integrate
        from sage.symbolic.integration.integral import definite_integral as integrate
        return integrate(expr, y, y_min, y_max) / (y_max - y_min)
    else:
        return expr

def compute_maxwell_stress_tensor_symbolic_ED(
    Ex, Ey, Ez,
    Dx, Dy, Dz,
    Hx, Hy, Hz,
    Bx, By, Bz,
):
    """
    Символьное вычисление компонент тензора Максвелла и их дивергенции

    Параметры:
    ----------
    Ex, Ey, Ez : символьные выражения компонент E поля
    Dx, Dy, Dz : символьные выражения компонент D поля
    Hx, Hy, Hz : символьные выражения компонент H поля
    Bx, By, Bz : символьные выражения компонент B поля

    Возвращает:
    -----------
    dict : словарь с символьными выражениями:
        - 'T_xx', 'T_xy', 'T_xz' : компоненты тензора
        - 'div_T_x' : дивергенция sum_k dT_xk/dx_k
        - dT_xx_dx = diff(T_xx, x)
        - dT_xy_dy = diff(T_xy, y)
        - dT_xz_dz = diff(T_xz, z)
    """

    from sage.functions.other import conjugate
    from sage.symbolic.constants import pi
    from sage.calculus.functional import diff


    # Скалярные произведения E·D и H·B
    ED = Ex*conjugate(Dx) + Ey*conjugate(Dy) + Ez*conjugate(Dz)
    # HB = Hx*conjugate(Bx) + Hy*conjugate(By) + Hz*conjugate(Bz)

    # Компоненты тензора Максвелла (формула 105.1 Тамма)
    # T_{ik} = 1/(8π) { E_i D_k + E_k D_i + H_i B_k + H_k B_i - δ_{ik}(E·D + H·B) }

    # Стандартный тензор Тамма (формула 105.1)
    # T_xx = 1/(8π) * { 2ExDx + 2HxBx - (ED + HB) }

    # T_xx (i=k=x, δ_xx = 1)
    T_xx = (1/(8*pi)) * (
        avg_over_y(2*Ex*conjugate(Dx) - (ED))/2
    )

    # T_xy (i=x, k=y, δ_xy = 0)
    T_xy = (1/(8*pi)) * (
        avg_over_y(Ex*conjugate(Dy) + Ey*conjugate(Dx))/2
    )

    # T_xz (i=x, k=z, δ_xz = 0)
    T_xz = (1/(8*pi)) * (
        avg_over_y(Ex*conjugate(Dz) + Ez*conjugate(Dx))/2
    )

    # Дивергенция: div T_x = ∂T_xx/∂x + ∂T_xy/∂y + ∂T_xz/∂z
    dT_xx_dx = diff(T_xx, x)
    dT_xy_dy = diff(T_xy, y)
    dT_xz_dz = diff(T_xz, z)

    div_T_x = dT_xx_dx + dT_xy_dy + dT_xz_dz

    return {
        'T_xx': T_xx,
        'T_xy': T_xy,
        'T_xz': T_xz,
        'div_T_x': div_T_x,
        'dT_xx_dx': dT_xx_dx,
        'dT_xy_dy': dT_xy_dy,
        'dT_xz_dz': dT_xz_dz,
    }

def compute_maxwell_stress_tensor_symbolic_HB(
    Ex, Ey, Ez,
    Dx, Dy, Dz,
    Hx, Hy, Hz,
    Bx, By, Bz,
):
    """
    Символьное вычисление компонент тензора Максвелла и их дивергенции

    Параметры:
    ----------
    Ex, Ey, Ez : символьные выражения компонент E поля
    Dx, Dy, Dz : символьные выражения компонент D поля
    Hx, Hy, Hz : символьные выражения компонент H поля
    Bx, By, Bz : символьные выражения компонент B поля

    Возвращает:
    -----------
    dict : словарь с символьными выражениями:
        - 'T_xx', 'T_xy', 'T_xz' : компоненты тензора
        - 'div_T_x' : дивергенция sum_k dT_xk/dx_k
        - dT_xx_dx = diff(T_xx, x)
        - dT_xy_dy = diff(T_xy, y)
        - dT_xz_dz = diff(T_xz, z)
    """

    from sage.functions.other import conjugate
    from sage.symbolic.constants import pi
    from sage.calculus.functional import diff


    # Скалярные произведения E·D и H·B
    # ED = Ex*conjugate(Dx) + Ey*conjugate(Dy) + Ez*conjugate(Dz)
    HB = Hx*conjugate(Bx) + Hy*conjugate(By) + Hz*conjugate(Bz)

    # Компоненты тензора Максвелла (формула 105.1 Тамма)
    # T_{ik} = 1/(8π) { E_i D_k + E_k D_i + H_i B_k + H_k B_i - δ_{ik}(E·D + H·B) }

    # Стандартный тензор Тамма (формула 105.1)
    # T_xx = 1/(8π) * { 2ExDx + 2HxBx - (ED + HB) }

    # T_xx (i=k=x, δ_xx = 1)
    T_xx = (1/(8*pi)) * (
        avg_over_y(2*Hx*conjugate(Bx) - (HB))/2
    )

    # T_xy (i=x, k=y, δ_xy = 0)
    T_xy = (1/(8*pi)) * (
        avg_over_y(Hx*conjugate(By) + Hy*conjugate(Bx))/2
    )

    # T_xz (i=x, k=z, δ_xz = 0)
    T_xz = (1/(8*pi)) * (
        avg_over_y(Hx*conjugate(Bz) + Hz*conjugate(Bx))/2
    )

    # Дивергенция: div T_x = ∂T_xx/∂x + ∂T_xy/∂y + ∂T_xz/∂z
    dT_xx_dx = diff(T_xx, x)
    dT_xy_dy = diff(T_xy, y)
    dT_xz_dz = diff(T_xz, z)

    div_T_x = dT_xx_dx + dT_xy_dy + dT_xz_dz

    return {
        'T_xx': T_xx,
        'T_xy': T_xy,
        'T_xz': T_xz,
        'div_T_x': div_T_x,
        'dT_xx_dx': dT_xx_dx,
        'dT_xy_dy': dT_xy_dy,
        'dT_xz_dz': dT_xz_dz,
    }

def compute_maxwell_stress_tensor_symbolic(
    Ex, Ey, Ez,
    Dx, Dy, Dz,
    Hx, Hy, Hz,
    Bx, By, Bz,
):
    """
    Символьное вычисление компонент тензора Максвелла и их дивергенции

    Параметры:
    ----------
    Ex, Ey, Ez : символьные выражения компонент E поля
    Dx, Dy, Dz : символьные выражения компонент D поля
    Hx, Hy, Hz : символьные выражения компонент H поля
    Bx, By, Bz : символьные выражения компонент B поля

    Возвращает:
    -----------
    dict : словарь с символьными выражениями:
        - 'T_xx', 'T_xy', 'T_xz' : компоненты тензора
        - 'div_T_x' : дивергенция sum_k dT_xk/dx_k
        - dT_xx_dx = diff(T_xx, x)
        - dT_xy_dy = diff(T_xy, y)
        - dT_xz_dz = diff(T_xz, z)
    """

    from sage.functions.other import conjugate
    from sage.symbolic.constants import pi
    from sage.calculus.functional import diff


    # Скалярные произведения E·D и H·B
    ED = Ex*conjugate(Dx) + Ey*conjugate(Dy) + Ez*conjugate(Dz)
    HB = Hx*conjugate(Bx) + Hy*conjugate(By) + Hz*conjugate(Bz)

    # Компоненты тензора Максвелла (формула 105.1 Тамма)
    # T_{ik} = 1/(8π) { E_i D_k + E_k D_i + H_i B_k + H_k B_i - δ_{ik}(E·D + H·B) }

    # Стандартный тензор Тамма (формула 105.1)
    # T_xx = 1/(8π) * { 2ExDx + 2HxBx - (ED + HB) }

    # T_xx (i=k=x, δ_xx = 1)
    T_xx = (1/(8*pi)) * (
        avg_over_y(2*Ex*conjugate(Dx) + 2*Hx*conjugate(Bx) - (ED + HB))/2
    )

    # T_xy (i=x, k=y, δ_xy = 0)
    T_xy = (1/(8*pi)) * (
        avg_over_y(Ex*conjugate(Dy) + Ey*conjugate(Dx) + Hx*conjugate(By) + Hy*conjugate(Bx))/2
    )

    # T_xz (i=x, k=z, δ_xz = 0)
    T_xz = (1/(8*pi)) * (
        avg_over_y(Ex*conjugate(Dz) + Ez*conjugate(Dx) + Hx*conjugate(Bz) + Hz*conjugate(Bx))/2
    )

    # Дивергенция: div T_x = ∂T_xx/∂x + ∂T_xy/∂y + ∂T_xz/∂z
    dT_xx_dx = diff(T_xx, x)
    dT_xy_dy = diff(T_xy, y)
    dT_xz_dz = diff(T_xz, z)

    div_T_x = dT_xx_dx + dT_xy_dy + dT_xz_dz

    return {
        'T_xx': T_xx,
        'T_xy': T_xy,
        'T_xz': T_xz,
        'div_T_x': div_T_x,
        'dT_xx_dx': dT_xx_dx,
        'dT_xy_dy': dT_xy_dy,
        'dT_xz_dz': dT_xz_dz,
    }
def compute_convective_stress_tensor_symbolic_PE(
    Ex, Ey, Ez,
    Dx, Dy, Dz,
    Hx, Hy, Hz,
    Bx, By, Bz,
):
    """
    Символьное вычисление компонент тензора Максвелла и их дивергенции

    Параметры:
    ----------
    Ex, Ey, Ez : символьные выражения компонент E поля
    Dx, Dy, Dz : символьные выражения компонент D поля
    Hx, Hy, Hz : символьные выражения компонент H поля
    Bx, By, Bz : символьные выражения компонент B поля

    Возвращает:
    -----------
    dict : словарь с символьными выражениями:
        - 'T_xx', 'T_xy', 'T_xz' : компоненты тензора
        - 'div_T_x' : дивергенция sum_k dT_xk/dx_k
        - dT_xx_dx = diff(T_xx, x)
        - dT_xy_dy = diff(T_xy, y)
        - dT_xz_dz = diff(T_xz, z)
    """

    from sage.functions.other import conjugate
    from sage.symbolic.constants import pi
    from sage.calculus.functional import diff


    # Скалярные произведения E·D и H·B
    EE = avg_over_y(Ex*conjugate(Ex) + Ey*conjugate(Ey) + Ez*conjugate(Ez))/2
    ED = avg_over_y(Ex*conjugate(Dx) + Ey*conjugate(Dy) + Ez*conjugate(Dz))/2
    HH = avg_over_y(Hx*conjugate(Hx) + Hy*conjugate(Hy) + Hz*conjugate(Hz))/2
    HB = avg_over_y(Hx*conjugate(Bx) + Hy*conjugate(By) + Hz*conjugate(Bz))/2

    # T_xx (i=k=x, δ_xx = 1)
    T_xx = (1/(8*pi)) * (
        ED - EE
    )

    # T_xy (i=x, k=y, δ_xy = 0)
    T_xy = (1/(8*pi)) * (
        0
    )

    # T_xz (i=x, k=z, δ_xz = 0)
    T_xz = (1/(8*pi)) * (
        0
    )

    # Дивергенция: div T_x = ∂T_xx/∂x + ∂T_xy/∂y + ∂T_xz/∂z
    dT_xx_dx = diff(T_xx, x)
    dT_xy_dy = diff(T_xy, y)
    dT_xz_dz = diff(T_xz, z)

    div_T_x = dT_xx_dx + dT_xy_dy + dT_xz_dz

    return {
        'T_xx': T_xx,
        'T_xy': T_xy,
        'T_xz': T_xz,
        'div_T_x': div_T_x,
        'dT_xx_dx': dT_xx_dx,
        'dT_xy_dy': dT_xy_dy,
        'dT_xz_dz': dT_xz_dz,
    }
def compute_convective_stress_tensor_symbolic_IH(
    Ex, Ey, Ez,
    Dx, Dy, Dz,
    Hx, Hy, Hz,
    Bx, By, Bz,
):
    """
    Символьное вычисление компонент тензора Максвелла и их дивергенции

    Параметры:
    ----------
    Ex, Ey, Ez : символьные выражения компонент E поля
    Dx, Dy, Dz : символьные выражения компонент D поля
    Hx, Hy, Hz : символьные выражения компонент H поля
    Bx, By, Bz : символьные выражения компонент B поля

    Возвращает:
    -----------
    dict : словарь с символьными выражениями:
        - 'T_xx', 'T_xy', 'T_xz' : компоненты тензора
        - 'div_T_x' : дивергенция sum_k dT_xk/dx_k
        - dT_xx_dx = diff(T_xx, x)
        - dT_xy_dy = diff(T_xy, y)
        - dT_xz_dz = diff(T_xz, z)
    """

    from sage.functions.other import conjugate
    from sage.symbolic.constants import pi
    from sage.calculus.functional import diff


    # Скалярные произведения E·D и H·B
    EE = avg_over_y(Ex*conjugate(Ex) + Ey*conjugate(Ey) + Ez*conjugate(Ez))/2
    ED = avg_over_y(Ex*conjugate(Dx) + Ey*conjugate(Dy) + Ez*conjugate(Dz))/2
    HH = avg_over_y(Hx*conjugate(Hx) + Hy*conjugate(Hy) + Hz*conjugate(Hz))/2
    HB = avg_over_y(Hx*conjugate(Bx) + Hy*conjugate(By) + Hz*conjugate(Bz))/2

    # T_xx (i=k=x, δ_xx = 1)
    T_xx = (1/(8*pi)) * (
        HB - HH
    )

    # T_xy (i=x, k=y, δ_xy = 0)
    T_xy = (1/(8*pi)) * (
        0
    )

    # T_xz (i=x, k=z, δ_xz = 0)
    T_xz = (1/(8*pi)) * (
        0
    )

    # Дивергенция: div T_x = ∂T_xx/∂x + ∂T_xy/∂y + ∂T_xz/∂z
    dT_xx_dx = diff(T_xx, x)
    dT_xy_dy = diff(T_xy, y)
    dT_xz_dz = diff(T_xz, z)

    div_T_x = dT_xx_dx + dT_xy_dy + dT_xz_dz

    return {
        'T_xx': T_xx,
        'T_xy': T_xy,
        'T_xz': T_xz,
        'div_T_x': div_T_x,
        'dT_xx_dx': dT_xx_dx,
        'dT_xy_dy': dT_xy_dy,
        'dT_xz_dz': dT_xz_dz,
    }
