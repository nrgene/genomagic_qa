#!/usr/bin/env python3
import sys
sys.path.append('/home/ariel/clients/my-tools/src/genomagic_utils')
import similarity_utils
import snp_markers_statistics

def compute_similarity_match_for_two_samples(sample1, sample2, raw_sim, snp_df , jar_file, sim_threshold, window_size,
                                             trim_len, informative):
    test_df = similarity_utils.create_dataframe_for_test(snp_df, sample1, sample2, informative)
    selected_samples_index = ((raw_sim['s1'] == sample1) & (raw_sim['s2'] == sample2) | (raw_sim['s2'] == sample1) & (
                raw_sim['s1'] == sample2))
    hap_sim = raw_sim.loc[selected_samples_index, ['chr', 'start', 'end']].astype(int).drop_duplicates()
    markers_sim = similarity_utils.compute_similarities_between_sampels_with_jar(test_df, jar_file, sim_threshold, window_size, trim_len)
    return compare_dataframes_of_similarities(hap_sim, markers_sim)


def compute_similarity_match_for_multiple_sampels(samples, raw_sim, snp_df , jar_file, sim_threshold, window_size, trim_len, max_major_allele_freq, min_samples_presence):
    hap_sim_unique_len = 0
    marker_sim_unique_len = 0
    common_sim_len = 0
    informative = snp_markers_statistics.compute_informative_markers(snp_df, max_major_allele_freq, min_samples_presence)
    n = len(samples)
    for i in range(n):
        sample1 = samples[i]
        print("comparing {}".format(sample1))
        for j in range(i+1, n):
            sample2 = samples[j]
            assert sample1 in snp_df.columns
            assert sample2 in snp_df.columns
            [a,b,c] = compute_similarity_match_for_two_samples(sample1, sample2, raw_sim, snp_df , jar_file,
                                                               sim_threshold, window_size, trim_len, informative)
            hap_sim_unique_len += a
            marker_sim_unique_len += b
            common_sim_len += c
    return [hap_sim_unique_len, marker_sim_unique_len, common_sim_len]































