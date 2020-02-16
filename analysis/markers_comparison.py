from scipy import signal
import numpy as np
import pandas as pd


# This function is performing one sided moving average of the data
# data - vector of 1 and 0
# win_len - window length
def get_sliding_window_average(data, win_len):
    x1 = signal.lfilter((1.0 / win_len) * np.ones(win_len), [1], data)
    return x1.round(2)


# This function gets indices vector, and split it where there are gaps
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


# This function gets seq - vector of 1/0 (match/mismatch)
# It returns the indices of the similarity regions
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


def get_matching_alleles(my_df, include_informative):
    n = my_df.shape[0]
    result = np.zeros((n,), dtype=bool)

    for allele in range(1, 3):
        informative_allele = (my_df['informative'] == allele) | (my_df['informative'] == 3) | (not include_informative)
        sample1_has_allele = (my_df['sample1'] == allele) | (my_df['sample1'] == 3)
        sample2_has_allele = (my_df['sample2'] == allele) | (my_df['sample2'] == 3)
        matchine_with_allele = informative_allele & sample1_has_allele & sample2_has_allele
        result[matchine_with_allele] = True
    return pd.Series(result)


def get_mismatching_alleles(my_df, include_informative):
    conflicting_markers = ((my_df['sample1'] == 1) & (my_df['sample2'] == 2)) | ((my_df['sample1'] == 2) & (my_df['sample2'] == 1))
    markers_are_informative = (my_df['informative'] > 0) | (not include_informative)
    return conflicting_markers & markers_are_informative


def get_match_and_mismatch_markers(my_df, include_informative):
    n = my_df.shape[0]
    match = get_matching_alleles(my_df, include_informative)
    mismatch = get_mismatching_alleles(my_df, include_informative)
    assert not (match & mismatch).any()
    match_array = np.zeros(n, dtype=np.int8)
    match_array[mismatch] = -1
    match_array[match] = 1
    return match_array


def get_similar_regions(v, min_p, min_len, trim_len):
    ind = np.where(v != 0)[0]
    b = v[ind]
    b[b == -1] = 0
    r = get_similar_regions_without_ignored_markers(b, min_p, min_len, trim_len)
    return [[ind[x[0]], ind[x[1]]] for x in r]


def compute_similarities_between_samples(my_df, min_p, min_len, trim_len):
    assert set(['sample1', 'sample2', 'informative', 'pos']).issubset(my_df.columns), my_df.columns
    match_array = get_match_and_mismatch_markers(my_df, True)
    match_pr = round((match_array == 1).mean()*100,1)
    mismatch_pr = round((match_array == -1).mean()*100,1)
    sim_reg_ind = get_similar_regions(match_array, min_p, min_len, trim_len)
    similarities = [[my_df['pos'][x[0]],my_df['pos'][x[1]]] for x in sim_reg_ind]
    return [similarities, match_pr, mismatch_pr]




#data = [0,0,1,1,1,1,0,0,1, 1, 1, 1, 0,0,1,1,1]
#x1 = get_sliding_window_average(data, 3)
#true_ind = np.where(x1 >= 0.9)[0]
#print(true_ind)
#true_regions = separate_regions_by_gap(true_ind)
#print(true_regions)
#x = separate_regions_by_gap(data)
#x = get_sliding_window_average(data, 3)
#print(x[2:6])
#print(x)