import os
import sys
import numpy as np
import pandas as pd
working_dir_path = os.getcwd()
genomagic_qa_repo_path = '/'.join(working_dir_path.split('/')[:-1])
sys.path.append(genomagic_qa_repo_path)
import vcf_utils.vcf_iterator as vi
vcf_metadata_columns_num = 9

def handle_sim_line(parts):
    haps_count = len(parts[4].split(','))-1
    samples_data = parts[vcf_metadata_columns_num:]
    samples_count  = len(samples_data)
    haps_indicator_mat = np.zeros((haps_count, samples_count), dtype=np.bool)
    hap_id_arr = []
    for i in range(samples_count):
        s = int(samples_data[i].split('|')[0])-1
        if s < haps_count:
            haps_indicator_mat[s,i] = True
    for i in range(haps_count):
        hap_id_arr.append('{}_{}_{}'.format(parts[0], parts[1], i))
    return [haps_indicator_mat, hap_id_arr]


def similarity_to_indicator_dataframe(sim_file):
    sim_ite = vi.VcfIterator(sim_file)
    headers = sim_ite.get_headers()
    samples = headers[9:]
    samples_num = len(samples)
    haps_mat = np.zeros((0, samples_num), dtype=np.bool)
    haps_id = []
    while sim_ite.has_next():
        l = sim_ite.get_next_line_splitted()
        [m, s] = handle_sim_line(l)
        haps_mat = np.concatenate((haps_mat, m), axis=0)
        haps_id = np.concatenate((haps_id, s), axis=0)
    print(haps_id)
    df = pd.DataFrame(haps_mat, columns=samples)
    df['hap_id'] = haps_id
    return df

#def print_indicator_data_frame(out_file, df):

sim_file = '/home/ariel/Downloads/hap_sim.vcf'
df = similarity_to_indicator_dataframe(sim_file)
print(df.T)

