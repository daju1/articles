from sage.plot.line import Line

def plot_r_of_sphere(plot_data, p, r_p_met, r_n_met, t_r_p, t_r_n):
    if r_p_met is True or r_n_met is True:
        (min_v, max_v) = get_min_max_of_data(plot_data)
    if r_p_met is True:
        L = Line([t_r_p, t_r_p], [min_v, max_v],{'alpha':1,'thickness':1,'rgbcolor':(1,0,0),'legend_label':''})
        p.add_primitive(L)
    if r_n_met is True:
        L = Line([t_r_n, t_r_n], [min_v, max_v],{'alpha':1,'thickness':1,'rgbcolor':(0,0,1),'legend_label':''})
        p.add_primitive(L)