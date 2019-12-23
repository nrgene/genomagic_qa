import sys
import numpy as np
genomagic_qa_repo_path = '/'.join(sys.argv[0].split('/')[:-2])
sys.path.append(genomagic_qa_repo_path)
import vcf_utils.vcf_iterator as vi
vcf_info_columns_num = 9


def vcf_line_to_mat(parts):
    samples_num = len(parts) - vcf_info_columns_num
    haps_num = parts[4].count(',')
    X = np.zeros((haps_num, samples_num), dtype=np.int8)
    Y = []
    for i in range(samples_num):
        curr_val = parts[vcf_info_columns_num + i]
        assert curr_val.count('|') == 1
        assert curr_val.count(':') == 2
        g = int(curr_val.split('|')[0])-1
        if g < haps_num:
            X[g,i] = 2
    for i in range(haps_num):
        Y.append('{}_{}_{}'.format(parts[0], parts[1],i))
    return X,Y


def sim_vcf_to_mat(hap_sim_file_name):
    sim_ite = vi.VcfIterator(hap_sim_file_name)
    a = sim_ite.get_headers()
    samples = a[vcf_info_columns_num:]
    samples_num = len(samples)
    hap_mat = np.zeros((0, samples_num), dtype=np.int8)
    hap_id = []
    while sim_ite.has_next():
        parts = sim_ite.get_next_line_splitted()
        [m, y] = vcf_line_to_mat(parts)
        hap_mat = np.concatenate((hap_mat, m), axis=0)
        hap_id = hap_id + y
    return [hap_mat, hap_id, samples]


def write_mat_to_csv(csv_out_name, hap_mat, hap_id, samples):
    total_haps_count = hap_mat.shape[0]
    samples_num = hap_mat.shape[1]
    assert samples_num == len(samples)
    assert total_haps_count == len(hap_id)
    csv_out = open(csv_out_name, 'w')
    csv_out.write('\"\"')
    for i in range(total_haps_count):
        csv_out.write(',\"{}\"'.format(hap_id[i]))
    csv_out.write(',\"Corrected.Strain\"\n')
    for i in range(samples_num):
        csv_out.write('\"{}\"'.format(i + 1))
        for j in range(total_haps_count):
            csv_out.write(',\"{}\"'.format(hap_mat[j, i]))
        csv_out.write(',\"{}\"\n'.format(samples[i]))
    csv_out.close()


def sim_vcf_to_csv(hap_sim_file_name, csv_out_name):
    [hap_mat, hap_id, samples] = sim_vcf_to_mat(hap_sim_file_name)
    write_mat_to_csv(csv_out_name, hap_mat, hap_id, samples)


hap_sim_file_name = '/mnt/ariel/genomagic/PSG-20/parents_progeny_merged.vcf'
csv_out_name = '/mnt/ariel/genomagic/PSG-20/parents_progeny_merged.vcf.csv'
sim_vcf_to_csv(hap_sim_file_name, csv_out_name)














