import sys
import numpy as np
genomagic_qa_repo_path = '/'.join(sys.argv[0].split('/')[:-2])
sys.path.append(genomagic_qa_repo_path)
import vcf_utils.vcf_iterator as vi
vcf_info_columns_num = 9


def vcf_line_to_mat(parts):
    samples_num = len(parts) - vcf_info_columns_num
    haps_num = parts[4].count(',') + 1
    X = np.zeros((haps_num, samples_num), dtype=np.int8)
    for i in range(samples_num):
        curr_val = parts[vcf_info_columns_num + i]
        assert curr_val.count('|') == 1
        assert curr_val.count(':') == 2
        g = int(curr_val.split('|')[0])-1
        assert g >= 0, curr_val
        if g >= 0 and g < haps_num:
            X[g,i] = 2
    curr_loc = '{}_{}'.format(parts[0],parts[1])
    #for i in range(haps_num-1):
    #    Y.append('{}_{}_HAP{}'.format(parts[0], parts[1],i+1))
    #Y.append('{}_{}_NAHAP'.format(parts[0], parts[1]))
    #non_zero_ind = X.sum(axis=1)>0
    #non_zero_haps = [Y[i] for i, val in enumerate(non_zero_ind) if val]
    return X, curr_loc


def fix_sample_names(vcf_sample_names):
    samples = [s.replace('ds12_', 'DS12-').replace('ds11_', 'DS11-') for s in vcf_sample_names]
    replace_samples = {'4j105_3_4': '4j105-3-4', "5m20_2_5_2":"5m20-2-5-2", "cl0j095_4_6":"cl0j095-4-6",
                       "cl0j173_6_8":"cl0j173-6-8", "hs6_3976": "hs6-3976", "ld00_3309":"ld00-3309",
                       "ld01_5907":"ld01-5907", "ld01_5907":"ld01-5907", "ld02_4485":"ld02-4485", "ld02_9050":"ld02-9050",
                       "lg00_3372": "lg00-3372", "lg03_2979":"lg03-2979", "lg03_3191":"lg03-3191",
                       "lg04_4717":"lg04-4717", "lg04_6000":"lg04-6000", "lg05_4292":"lg05-4292",
                       "lg05_4317":"lg05-4317", "lg05_4464":"lg05-4464", "lg05_4832":"lg05-4832",
                       "lg90_2550":"lg90-2550", "lg92_1255":"lg90-2550", "lg92_1255":"lg92-1255",
                       "lg94_1128":"lg94-1128", "lg94_1906":"lg94-1906", "lg97_7012":"lg97-7012",
                       "lg98_1605":"lg98-1605", "parent_ia3023":"parent-ia3023", "s06_13640":"s06-13640",
                       "tn05_3027":"tn05-3027", "u03_100612":"u03-100612"}
    fixed_sampels = [replace_samples[x] if x in replace_samples else x for x in samples]
    return fixed_sampels


def sim_vcf_to_mat(hap_sim_file_name, max_haps_num):
    sim_ite = vi.VcfIterator(hap_sim_file_name)
    a = sim_ite.get_headers()
    #samples = a[vcf_info_columns_num:]
    samples = fix_sample_names(a[vcf_info_columns_num:]) #[s.replace('ds12_', 'DS12-').replace('ds11_', 'DS11-') for s in a[vcf_info_columns_num:]]
    samples_num = len(samples)
    hap_mat = np.zeros((max_haps_num, samples_num), dtype=np.int8)
    counter = 0
    hap_id = []
    while sim_ite.has_next():
        parts = sim_ite.get_next_line_splitted()
        m, curr_loc = vcf_line_to_mat(parts)
        row_num = m.shape[0]
        assert m.shape[1] == samples_num
        assert counter+row_num < max_haps_num
        hap_mat[counter: counter+row_num,:] = m
        hap_id.append(curr_loc)
        counter += m.shape[0]
    return [hap_mat[0:counter, :], hap_id, samples]


def write_mat_to_csv(csv_out_name, hap_mat, hap_id, samples):
    total_haps_count = hap_mat.shape[0]
    samples_num = hap_mat.shape[1]
    assert samples_num == len(samples)
    assert total_haps_count >= len(hap_id)
    csv_out = open(csv_out_name, 'w')
    csv_out.write('\"\"')
    vcf_lines_num = len(hap_id)
    for i in range(vcf_lines_num):
        csv_out.write(',\"{}\"'.format(hap_id[i]))
    for i in range(vcf_lines_num, total_haps_count):
        csv_out.write(',\"NA\"')
    csv_out.write(',\"Corrected.Strain\"\n')
    for i in range(samples_num):
        csv_out.write('\"{}\"'.format(i + 1))
        for j in range(total_haps_count):
            csv_out.write(',{}'.format(hap_mat[j, i]))
        csv_out.write(',\"{}\"\n'.format(samples[i]))
    csv_out.close()


def sim_vcf_to_csv(hap_sim_file_name, csv_out_name):
    [hap_mat, hap_id, samples] = sim_vcf_to_mat(hap_sim_file_name)
    write_mat_to_csv(csv_out_name, hap_mat, hap_id, samples)



#hap_sim_file_name = sys.argv[1]#'/mnt/ariel/genomagic/PSG-20/parents_progeny_merged.vcf'
#csv_out_name = sys.argv[2]#'/mnt/ariel/genomagic/PSG-20/parents_progeny_merged.vcf.csv'
#sim_vcf_to_csv(hap_sim_file_name, csv_out_name)













