def get_min_max_of_data(data):
    min = sys.float_info.max;
    max = -sys.float_info.max;
    for (t, v) in data:
        if min > v:
            min = v
        if max < v:
            max = v
    return (min, max)