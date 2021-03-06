import sys
import os
import pandas as pd
working_dir_path = os.getcwd()
genomagic_qa_repo_path = '/'.join(working_dir_path.split('/')[:-1])
sys.path.append(genomagic_qa_repo_path)


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


def get_subset_of_similarities(hap_sim1, sample1, sample2):
    I = ((hap_sim1['s1'] == sample1) & (hap_sim1['s2'] == sample2)) | ((hap_sim1['s2'] == sample1) & (hap_sim1['s1'] == sample2))
    if I.sum() == 0:
        return pd.DataFrame(columns=['chr', 'start', 'end'])
    else:
        df1 = hap_sim1.loc[I,['chr', 'start', 'end']]
        return df1.astype(int).drop_duplicates()


def get_similarities_comparison_by_samples_pair(hap_sim1, hap_sim2, samples):
    assert set(['s1', 's2', 'chr', 'start', 'end']).issubset(hap_sim1.columns), hap_sim1.columns
    assert set(['s1', 's2', 'chr', 'start', 'end']).issubset(hap_sim2.columns), hap_sim2.columns
    arr = []
    samples_num = len(samples)
    for i in range(samples_num):
        for j in range(i + 1, samples_num):
            sample1 = samples[i]
            sample2 = samples[j]
            df1 = get_subset_of_similarities(hap_sim1, sample1, sample2)
            df2 = get_subset_of_similarities(hap_sim2, sample1, sample2)
            [only_in_1_len, only_in_2_len, shared_length] = compare_dataframes_of_similarities(
                df1, df2)
            arr.append([sample1, sample2, only_in_1_len, only_in_2_len, shared_length])
    return pd.DataFrame(arr, columns=['sample1', 'sample2', 'unique in 1', 'unique in 2', 'both'])


def compute_similarity_match_for_multiple_sampels_in_df(hap_sim1, hap_sim2, samples):
    df = get_similarities_comparison_by_samples_pair(hap_sim1, hap_sim2, samples)
    hap_sim1_uniq_len = df['unique in 1'].sum()
    hap_sim2_uniq_len = df['unique in 2'].sum()
    shared_len = df['both'].sum()
    return [hap_sim1_uniq_len, hap_sim2_uniq_len, shared_len]








#[fn1, fp1, intersect1] = compute_similarity_match_for_multiple_sampels_in_df(hap_sim1, hap_sim2, samples)
#print("{}, {}, {}".format(fn1/1000000, fp1/1000000, intersect1/1000000))

#b73
#cml247
#mo17
#w22
#ki3


