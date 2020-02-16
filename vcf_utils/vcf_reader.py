import numpy as np
import pandas as pd
vcf_info_columns_num = 9


def parse_snp_vcf_line(parts, my_dict):
    my_chr = int(parts[0])
    my_pos = int(parts[1])
    my_id = parts[2]
    samples_num = len(parts) - vcf_info_columns_num
    my_arr = np.zeros((samples_num,), dtype=np.int8)
    for i in range(samples_num):
        curr_val = parts[i + vcf_info_columns_num].split(':')[0]
        assert curr_val in my_dict
        my_arr[i] = my_dict[curr_val]
    return [my_chr, my_pos, my_id, my_arr]


def read_snp_vcf_from_list(vcf_lines_list, max_haps_num):
    my_dict = {'1|1': 1, '2|2': 2, '1|2': 3, '2|1': 3, '3|3': 0, '4|4': 0}
    curr_vcf_line = vcf_lines_list.pop(0)
    assert curr_vcf_line[0] == '#', curr_vcf_line[0]
    while curr_vcf_line[0] == '#':
        prev = curr_vcf_line
        curr_vcf_line = vcf_lines_list.pop(0)
    headers = prev.split('\t')
    samples_list = headers[vcf_info_columns_num:]
    samples_num = len(samples_list)
    hap_mat = np.zeros((max_haps_num, samples_num), dtype=np.int8)
    my_chr = []
    my_pos = []
    my_id = []
    counter = 0
    while vcf_lines_list:
        assert counter < max_haps_num
        [c,p,i,a] = parse_snp_vcf_line(curr_vcf_line.split('\t'), my_dict)
        hap_mat[counter,:] = a
        my_chr.append(c)
        my_pos.append(p)
        my_id.append(i)
        curr_vcf_line = vcf_lines_list.pop(0)
        counter += 1
    data = {'chr': my_chr, 'pos': my_pos, 'id':my_id}
    for i in range(samples_num):
        data[samples_list[i]] = hap_mat[:counter,i]
    return pd.DataFrame(data)


def read_hap_vcf_from_list(vcf_lines_list, max_haps_num):
    my_dict = {'1|1': False, '2|2': True}
    curr_vcf_line = vcf_lines_list.pop(0)
    assert curr_vcf_line[0] == '#', curr_vcf_line[0]
    while curr_vcf_line[0] == '#':
        prev = curr_vcf_line
        curr_vcf_line = vcf_lines_list.pop(0)
    headers = prev.split('\t')
    samples_list = headers[vcf_info_columns_num:]
    samples_num = len(samples_list)
    hap_mat = np.zeros((max_haps_num, samples_num), dtype=bool)
    my_chr = []
    my_pos = []
    my_id = []
    counter = 0
    while vcf_lines_list:
        assert counter < max_haps_num
        [c,p,i,a] = parse_snp_vcf_line(curr_vcf_line.split('\t'), my_dict)
        hap_mat[counter,:] = a
        my_chr.append(c)
        my_pos.append(p)
        my_id.append(i)
        curr_vcf_line = vcf_lines_list.pop(0)
        counter += 1
    data = {'chr': my_chr, 'pos': my_pos, 'id':my_id}
    for i in range(samples_num):
        data[samples_list[i]] = hap_mat[:counter,i]
    return pd.DataFrame(data)

