from scipy import signal
import numpy as np


# This function is performing one sided moving average of the data
def get_sliding_window_average(data, win_len):
    x1 = signal.lfilter((1.0 / win_len) * np.ones(win_len), [1], data)
    return x1.round(2)


# This function is
def separate_regions_by_gap(v):
    res = []
    n = len(v)
    if n == 0:
        return []
    g = np.where(np.diff(v) > 1)[0]
    g = np.append(g, n-1)
    res.append([0, g[0]])
    if g.size > 1:
        for e in g[1:]:
            vec = [res[-1][1]+1, e]
            res.append(vec)
    return [[v[x[0]], v[x[1]]] for x in res]


def get_similar_regions_without_ignored_markers(seq, min_p, min_len, trim_len):
    x1 = get_sliding_window_average(seq, min_len)
    true_ind = np.where(x1 >= min_p)[0]
    true_regions = separate_regions_by_gap(true_ind)
    my_regions = []
    for t in true_regions:
        start_point = max(t[0]-min_len+1 + trim_len,0)
        end_point = t[1]-trim_len
        my_regions.append([start_point, end_point])
    return my_regions



data = [0,0,1,1,1,1,0,0,1, 1, 1, 1, 0,0,1,1,1]
x1 = get_sliding_window_average(data, 3)
true_ind = np.where(x1 >= 0.9)[0]
print(true_ind)
true_regions = separate_regions_by_gap(true_ind)
print(true_regions)
#x = separate_regions_by_gap(data)
#x = get_sliding_window_average(data, 3)
#print(x[2:6])
#print(x)