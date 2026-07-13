def calc_extrap_len(kz_range, sz_range, nk, ns):
    import numpy as np
    extrap_len = np.sqrt(
        ((kz_range[1] - kz_range[0])/nk)**2 +
        ((sz_range[1] - sz_range[0])/ns)**2)
    return extrap_len

# deprecated
def init_det_plot2d(kz_range, sz_range, nk, ns):
    import numpy as np

    kz_min, kz_max = kz_range
    sz_min, sz_max = sz_range

    kz_linspace = np.linspace(real128(kz_min), real128(kz_max), nk)
    sz_linspace = np.linspace(real128(sz_min), real128(sz_max), ns)

    kz_list = kz_linspace.tolist()
    sz_list = sz_linspace.tolist()

    kz_grid, sz_grid = np.meshgrid(kz_linspace, sz_linspace)

    u = kz_grid * np.nan
    v = sz_grid * np.nan

    shl = kz_grid * np.nan
    sel = sz_grid * np.nan
    shr = kz_grid * np.nan
    ser = sz_grid * np.nan

    return kz_grid, sz_grid, kz_list, sz_list, u, v, shl, sel, shr, ser
