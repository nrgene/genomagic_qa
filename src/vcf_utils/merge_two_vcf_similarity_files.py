import sys
import vcf_pair_iterator
vcf_info_columns_num = 9

def get_index_map(samples1, samples2):
    n1 = len(samples1)
    index_from_progeny_to_parents = {}
    for i in range(n1):
        curr_sample = samples1[i]
        if curr_sample in samples2:
            index_from_progeny_to_parents[str(i+1)] = samples2.index(curr_sample)
    index_from_progeny_to_parents[str(len(index_from_progeny_to_parents)+1)] = -1
    return index_from_progeny_to_parents


def get_curr_value(index_in_parents_vcf, parts2, nahap_ind, seg_len):
    if index_in_parents_vcf == -1:
        return '{}|{}:0,0,0:{}'.format(nahap_ind, nahap_ind, seg_len)
    else:
        new_val = parts2[vcf_info_columns_num + index_in_parents_vcf].split(':')
        assert len(new_val) == 3
        return '{}:{}:{}'.format(new_val[0], new_val[1], seg_len)


def get_line_info(chr_name, start, end, haps_list):
    curr_str = '{}\t{}\t.\tN\t{}\t1\tPASS\tHS;END={}\tGT:CO1:LN'.format(chr_name, start, haps_list, end)
    return curr_str


def get_line_data(prev_data, index_from_progeny_to_parents, target_data_line, seg_len):
    curr_str=''
    nahap_ind = target_data_line[4].count(',')
    for g in prev_data:
        index_in_parents_vcf = index_from_progeny_to_parents[g.split('|')[0]]
        tmp_str = get_curr_value(index_in_parents_vcf, target_data_line, nahap_ind, seg_len)
        curr_str = '{}\t{}'.format(curr_str, tmp_str)
    return curr_str


def my_simple_func(line1, line2, index_from_progeny_to_parents):
    [chr1, start1, end1] = vcf_pair_iterator.parse_vcf_line(line1)
    [chr2, start2, end2] = vcf_pair_iterator.parse_vcf_line(line1)
    assert chr1 == chr2
    assert start1 < end2
    assert start2 < end1
    my_start = max(start1, start2)
    my_end = min(end1, end2)
    parts1 = line1.split('\t')
    parts2 = line2.split('\t')
    info_part = get_line_info(chr1, my_start, my_end, parts2[4])
    data_part = get_line_data(parts1[vcf_info_columns_num:], index_from_progeny_to_parents, parts2, my_end - my_start + 1)
    print('{}{}'.format(info_part, data_part))


progeny_vcf_file = sys.argv[1]
parents_vcf_file = sys.argv[2]
first_file = open(progeny_vcf_file, 'r')
second_file = open(parents_vcf_file, 'r')
[line1, header1] = vcf_pair_iterator.get_first_non_header_line(first_file)
[line2, header2] = vcf_pair_iterator.get_first_non_header_line(second_file)
samples1 = header1.split('\t')[vcf_info_columns_num:]
samples2 = header2.split('\t')[vcf_info_columns_num:]
index_from_progeny_to_parents = get_index_map(samples1, samples2)
vcf_pair_iterator.iterate_two_vcf(first_file, second_file, line1, line2, my_simple_func, index_from_progeny_to_parents)
first_file.close()
second_file.close()


#
