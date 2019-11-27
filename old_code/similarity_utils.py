import pandas as pd
import numpy as np
from scipy import signal
import sys
import os
from os import path
from subprocess import Popen, PIPE, STDOUT
sys.path.append('../vcf_utils')
import snp_markers_statistics


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


def get_sliding_window_average(data, win_len):
    x1 = signal.lfilter((1.0 / win_len) * np.ones(win_len), [1], data)
    return x1.round(2)


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


def create_dataframe_for_test(df, sample1, sample2, informative):
    df_test = pd.DataFrame(
        {'chr':df['chr'], 'pos': df['pos'], 'informative': informative, 'sample1': df[sample1], 'sample2': df[sample2]})
    return df_test


def get_match_and_mismatch_markers(my_df, include_informative):
    n = my_df.shape[0]
    match = get_matching_alleles(my_df, include_informative)
    mismatch = get_mismatching_alleles(my_df, include_informative)
    assert not (match & mismatch).any()
    match_array = np.zeros(n, dtype=np.int8)
    match_array[mismatch] = -1
    match_array[match] = 1
    return match_array


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


def get_probability_of_random_similarity_sequence(total_markers_num, p, win_len, min_th, trials_num):
    success_num = 0.0
    for i in range(trials_num):
        x = 1.0 * (np.random.rand(total_markers_num, ) > p)
        y = get_sliding_window_average(x, win_len)
        if any(y >= min_th):
            success_num += 1
    return success_num / trials_num


def get_rand_markers_with_same_probability(count1, count2):
    p1 = count1 / (count1 + count2)
    r = np.random.rand(p1.shape[0],)
    return (r > p1) + 1


def create_simulated_data(df, samples_num, max_major_allele_freq, min_samples_presence):
    informative = snp_markers_statistics.compute_informative_markers(df, max_major_allele_freq, min_samples_presence)
    [count1, count2] = snp_markers_statistics.get_alleles_count(df)
    df_sim = pd.DataFrame({'pos': df['pos'], 'informative': informative})
    for i in range(samples_num):
        v = get_rand_markers_with_same_probability(count1, count2)
        sample_name = "sample_{}".format(i)
        df_sim[sample_name] = v
    return df_sim


def get_similarities_from_two_sampels(df, sample1, sample2, informative, min_p, min_len, trim_len):
    df_temp = create_dataframe_for_test(df, sample1, sample2, informative)
    [arr, pos, neg] = compute_similarities_between_samples(df_temp, min_p, min_len, trim_len)
    return [sum([x[1] - x[0] for x in arr]), pos, neg]


def compute_mean_similarity_from_simulation(df, simulations_num, max_major_allele_freq, min_samples_presence, min_p, min_len, trim_len):
    sim_data = create_simulated_data(df, simulations_num, max_major_allele_freq, min_samples_presence)
    comparisons_num = 0
    total_sim_len = 0.0
    total_match = 0.0
    total_mismatch = 0.0
    for i in range(simulations_num):
        for j in range(i + 1, simulations_num):
            comparisons_num += 1
            [curr_sim_len, pos, neg] = get_similarities_from_two_sampels(sim_data, "sample_{}".format(i), "sample_{}".format(j),
                                                               sim_data['informative'], min_p, min_len, trim_len)
            total_sim_len += curr_sim_len
            total_match += pos
            total_mismatch += neg
    mean_sim_len = total_sim_len / comparisons_num
    mean_match = total_match / comparisons_num
    mean_mismatch = total_mismatch / comparisons_num
    return [mean_sim_len, mean_match, mean_mismatch]


def percentage_of_matching_markers(df, sample1, sample2, informative):
    test_df = create_dataframe_for_test(df, sample1, sample2, informative)
    match_vec = get_match_and_mismatch_markers(test_df, True)
    match_mean = (match_vec == 1).mean()
    mismatch_mean = (match_vec == -1).mean()
    pr = match_mean/(match_mean + mismatch_mean)
    return pr


def run_similarity_java(jar_file, input_file, similarityThreshold, slidingWindowSize, trimming):
    assert path.exists(jar_file)
    cmd = 'java -jar {} SlidingWindow --inputFile {} --similarityThreshold {} --slidingWindowSize {} --trimming {}'\
        .format(jar_file, input_file, similarityThreshold, slidingWindowSize, trimming)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
    output = str(p.stdout.read()).replace("\\t", "\t").replace("\\n", "\n")
    output_lines = output[2:-2].split("\n")[1:]
    data = [x.split(',') for x in output_lines]
    if len(data) == 0:
        return pd.DataFrame(columns=['chr', 'start', 'end'])
    assert len(data) > 0
    df = pd.DataFrame(data)
    assert df.shape[1] == 5, df.shape
    df.columns = ["chr", "start", "end", "score", "ibd"]
    df = df.loc[df["ibd"] == "true", ['chr', 'start', 'end']].astype(int).drop_duplicates()
    return df


def write_temp_file_of_markers_comparison(filename, my_df, match_array):
    assert my_df.shape[0] == len(match_array)
    assert 'chr' in my_df.columns
    assert 'pos' in my_df.columns
    n = my_df.shape[0]
    print(n)
    f = open(filename, "w")
    for i in range(n):
        if match_array[i] != 0:
            my_chr = my_df['chr'][i]
            my_pos = my_df['pos'][i]
            my_val = 'true' if match_array[i] == 1 else 'false'
            f.write("{},{},{},{},0\n".format(my_chr, my_pos, my_pos + 1, my_val))
    f.close()


def compute_similarities_between_sampels_with_jar(my_df, jar_file, similarityThreshold, slidingWindowSize, trimming):
    assert set(['sample1', 'sample2', 'informative', 'pos']).issubset(my_df.columns), my_df.columns
    match_array = get_match_and_mismatch_markers(my_df, True)
    my_temp_file = 'temp_file_for_markers'
    write_temp_file_of_markers_comparison(my_temp_file, my_df, match_array)
    true_similarities = run_similarity_java(jar_file, my_temp_file, similarityThreshold, slidingWindowSize, trimming)
    #os.remove(my_temp_file)
    return true_similarities









