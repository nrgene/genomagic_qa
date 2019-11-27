#!/usr/bin/env python3
import sys
sys.path.append('/home/ariel/clients/my-tools/src/genomagic_utils')
import similarity_utils
import snp_markers_statistics

def advance_in_list(my_list, prev_element_end):
    if len(my_list) > 0:
        seg = my_list.pop(0)
        assert seg[0] < seg[1], "expecting seg[0]={} < seg[1]={}".format(seg[0] , seg[1])
        assert seg[0] > prev_element_end, "expecting seg[0]={} > prev_element_end={}".format(seg[0] , prev_element_end)
        return [seg[0], seg[1]]
    else:
        return [None, None]


def handle_partially_overlap(first_start, first_end, second_start, second_end, list_first):
    assert first_start <= second_start
    assert second_start <= first_end
    assert first_end <= second_end
    overlap_len = first_end - second_start
    exists_only_in_first_len = second_start - first_start
    second_start = first_end
    [first_start, first_end] = advance_in_list(list_first, first_end)
    return [first_start, first_end, second_start, exists_only_in_first_len, overlap_len]


# if a contains b, then advance second, and set a_start to be the end of b
def handle_when_a_contains_b(a_start, a_end, b_start, b_end, list_b):
    assert a_start <= b_start
    assert a_end >= b_end, "a ({}-{}) contains b({}-{}) ??".format(a_start, a_end, b_start, b_end)    
    overlap_len = b_end - b_start
    unique_to_a = b_start - a_start
    assert unique_to_a >= 0
    a_start = b_end
    [b_start, b_end] = advance_in_list(list_b, b_end)
    return [a_start, b_start, b_end, unique_to_a, overlap_len]


# if first  overlaps with second
def compare_lists_inrtersect(list1, list2):
    only_in_list1_length = 0
    only_in_list2_length = 0
    lists_overlap_length = 0
    [start1, end1] = advance_in_list(list1, -1)
    [start2, end2] = advance_in_list(list2, -1)
    while (start1 is not None) or (start2 is not None):

        # if 1 is none, advance in 2
        if start1 is None:
            curr_len = end2 - start2
            assert curr_len >= 0
            only_in_list2_length += curr_len
            [start2, end2] = advance_in_list(list2, end2)

        # if 2 is none, advance in 1
        elif start2 is None:
            curr_len = end1 - start1
            assert curr_len >= 0
            only_in_list1_length += curr_len
            [start1, end1] = advance_in_list(list1, end1)

        # if 2 is before 1, no overlap
        elif start1 >= end2:
            curr_len = end2 - start2
            assert curr_len >= 0
            only_in_list2_length += curr_len
            [start2, end2] = advance_in_list(list2, end2)

        # if 1 is before 2, no overlap
        elif start2 >= end1:
            curr_len = end1 - start1
            assert curr_len >= 0
            only_in_list1_length += curr_len
            [start1, end1] = advance_in_list(list1, end1)

        # if seg1 before seg2 with overlap, advance seg1, and update start of seg2
        elif start1 <= start2 and end1 <= end2:
            [start1, end1, start2, unique_len, overlap_len] = handle_partially_overlap(start1, end1, start2, end2, list1)
            assert unique_len >= 0
            assert overlap_len >= 0
            only_in_list1_length += unique_len
            lists_overlap_length += overlap_len

        # if seg2 before seg1 with overlap, advance seg2, and update start of seg1
        elif start1 > start2 and end1 > end2:
            [start2, end2, start1, unique_len, overlap_len] = handle_partially_overlap(start2, end2, start1, end1, list2)
            assert unique_len >= 0
            assert overlap_len >= 0
            only_in_list2_length += unique_len
            lists_overlap_length += overlap_len

        # seg1 contains seg2
        elif start1 <= start2 and end1 > end2:
            [start1, start2, end2, unique_len, overlap_len] = handle_when_a_contains_b(start1, end1, start2, end2,
                                                                                       list2)
            assert unique_len >= 0
            assert overlap_len >= 0
            only_in_list1_length += unique_len
            lists_overlap_length += overlap_len

        # seg2 contains seg1
        else:
            [start2, start1, end1, unique_len, overlap_len] = handle_when_a_contains_b(start2, end2, start1, end1,
                                                                                       list1)
            assert unique_len >= 0
            assert overlap_len >= 0
            only_in_list2_length += unique_len
            lists_overlap_length += overlap_len
    assert only_in_list1_length >= 0, only_in_list1_length
    assert only_in_list2_length >= 0, only_in_list2_length
    assert lists_overlap_length >= 0, lists_overlap_length

    return [only_in_list1_length, only_in_list2_length, lists_overlap_length]


def merge_list_elements(my_list):
    if len(my_list) == 0:
        return my_list
    merged_list = []
    curr_val = [my_list[0][0], my_list[0][1]]
    n = len(my_list)
    for i in range(1,n):
        assert my_list[i][0] <= my_list[i][1],"expecting my_list[i][0] =  {} < my_list[i][1] = {}".format(my_list[i][0] , my_list[i][1])
        #assert my_list[i][0] >= my_list[i-1][1], "expecting my_list[i] = {} >= my_list[i-1] = {}".format(my_list[i] , my_list[i-1])
        if my_list[i][0] <= my_list[i-1][1]:
            curr_val[1] = max(my_list[i][1],my_list[i][1])
        else:
            merged_list.append(curr_val)
            curr_val = [my_list[i][0], my_list[i][1]]
    merged_list.append(curr_val)
    return merged_list


def compare_dataframes_of_similarities(df1, df2):
    assert set(['chr','start', 'end']) == set(df1.columns), df1.columns
    assert set(['chr','start', 'end']) == set(df2.columns), df2.columns
    all_chr = df1['chr'].append(df2['chr']).unique()
    only_in_1_len = 0
    only_in_2_len = 0
    shared_length = 0
    for c in all_chr:
        list1 = merge_list_elements(list(df1.loc[df1["chr"] == c, ['start', 'end']].values))
        list2 = merge_list_elements(list(df2.loc[df2["chr"] == c, ['start', 'end']].values))
        [only_in_list1_length, only_in_list2_length, lists_overlap_length] = compare_lists_inrtersect(list1, list2)
        only_in_1_len += only_in_list1_length
        only_in_2_len += only_in_list2_length
        shared_length += lists_overlap_length
    return [only_in_1_len, only_in_2_len, shared_length]


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


def get_subset_of_similarities(hap_sim1, sample1, sample2):
    I = ((hap_sim1['s1'] == sample1) & (hap_sim1['s2'] == sample2)) | ((hap_sim1['s2'] == sample1) & (hap_sim1['s1'] == sample2))
    df1 = hap_sim1.loc[I,['chr', 'start', 'end']]
    return df1.astype(int).drop_duplicates()


def compute_similarity_match_for_multiple_sampels_in_df(hap_sim1, hap_sim2, samples):
    samples_num = len(samples)
    hap_sim1_uniq_len = 0
    hap_sim2_uniq_len = 0
    shared_len = 0
    for i in range(samples_num):
        for j in range(i + 1, samples_num):
            sample1 = samples[i]
            sample2 = samples[j]
            df1 = get_subset_of_similarities(hap_sim1, sample1, sample2)
            df2 = get_subset_of_similarities(hap_sim2, sample1, sample2)
            [only_in_1_len, only_in_2_len, shared_length] = compare_dataframes_of_similarities(
                df1, df2)
            hap_sim1_uniq_len += only_in_1_len
            hap_sim2_uniq_len += only_in_2_len
            shared_len += shared_length
    return [hap_sim1_uniq_len, hap_sim2_uniq_len, shared_len]































