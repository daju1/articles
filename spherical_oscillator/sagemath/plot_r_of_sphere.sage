from sage.plot.line import Line

def plot_r_of_sphere(plot_data, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max):
    return

    (min_v, max_v) = get_min_max_of_data(plot_data)
    if r_p_min_met is True:
        L = Line([t_r_p_min, t_r_p_min], [min_v, max_v],{'alpha':1,'thickness':1,'rgbcolor':(1,0,0),'legend_label':''})
        p.add_primitive(L)
    if r_p_max_met is True:
        L = Line([t_r_p_max, t_r_p_max], [min_v, max_v],{'alpha':1,'thickness':1,'rgbcolor':(1,0,0),'legend_label':''})
        p.add_primitive(L)
    if r_n_min_met is True:
        L = Line([t_r_n_min, t_r_n_min], [min_v, max_v],{'alpha':1,'thickness':1,'rgbcolor':(0,0,1),'legend_label':''})
        p.add_primitive(L)
    if r_n_max_met is True:
        L = Line([t_r_n_max, t_r_n_max], [min_v, max_v],{'alpha':1,'thickness':1,'rgbcolor':(0,0,1),'legend_label':''})
        p.add_primitive(L)