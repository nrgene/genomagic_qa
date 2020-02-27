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


# return vector that is iff marker has selected_allele in both samples and is informative
def get_match_by_specific_allele(sample1_markers, sample2_markers, marker_informative, selected_allele):
    both_alleles = 3
    allele_is_informative = (marker_informative == selected_allele) | (marker_informative == both_alleles)
    sample1_has_allele = (sample1_markers == selected_allele) | (sample1_markers == both_alleles)
    sample2_has_allele = (sample2_markers == selected_allele) | (sample2_markers == both_alleles)
    return allele_is_informative & sample1_has_allele & sample2_has_allele


# return vector that is iff both samples have the same allele and its informative
def get_matching_alleles(my_df, include_informative):
    both_alleles = 3
    marker_informative = my_df['informative']
    if not include_informative:
        ind_not_full_informative = marker_informative != both_alleles
        marker_informative[ind_not_full_informative] = both_alleles
    match_by_allele1 = get_match_by_specific_allele(my_df['sample1'], my_df['sample2'], marker_informative, 1)
    match_by_allele2 = get_match_by_specific_allele(my_df['sample1'], my_df['sample2'], marker_informative, 2)
    matching_alleles = match_by_allele1 | match_by_allele2
    return matching_alleles

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


def compute_informative_markers(snp_markers, max_major_allele_freq, min_samples_presence):
    df = snp_markers.loc[:, snp_markers.columns[3:]]
    count1 = (df == 1).sum(axis=1)
    count2 = (df == 2).sum(axis=1)
    count_all = count1 + count2
    informative = (1 * ((count1 / count_all) < max_major_allele_freq)
                   + 2 * ((count2 / count_all) < max_major_allele_freq)) \
                  * (count_all > min_samples_presence)
    return informative


def get_markers_of_two_samples_for_comparison(df, sample1, sample2, informative):
    df_test = pd.DataFrame(
        {'chr':df['chr'], 'pos': df['pos'], 'informative': informative, 'sample1': df[sample1], 'sample2': df[sample2]})
    return df_test


def add_snp_similarities_to_dataframe(snp_similarities, snp_markers, chr_id, informative, sample1, sample2, min_p, win_len, trim_len):
    pair_markers = get_markers_of_two_samples_for_comparison(snp_markers, sample1, sample2, informative)
    [snp_sim_pair, match_pr, mismatch_pr] = compute_similarities_between_samples(pair_markers, min_p, win_len, trim_len)
    for s in snp_sim_pair:
        snp_similarities = snp_similarities.append({'s1': sample1, 's2': sample2, 'chr':chr_id, 'start': s[0], 'end':s[1]}, ignore_index=True)
    return snp_similarities


def get_all_pairwise_similarities_snp(snp_markers, samples_list, chr_id, min_p, win_len, trim_len , max_major_allele_freq, min_samples_presence):
    informative = compute_informative_markers(snp_markers, max_major_allele_freq, min_samples_presence)
    samples_count = len(samples_list)
    snp_similarities = pd.DataFrame(columns=['s1', 's2', 'chr','start', 'end'])
    for i in range(samples_count):
        #print("i={}".format(i))
        for j in range(i + 1, samples_count):
            sample1 = samples_list[i]
            sample2 = samples_list[j]
            snp_similarities = add_snp_similarities_to_dataframe(snp_similarities, snp_markers,chr_id,  informative, sample1,
                                                                 sample2, min_p, win_len, trim_len)
    return snp_similarities

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