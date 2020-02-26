import os
import sys
genomagic_abs_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append('{}/vcf_utils'.format(genomagic_abs_path))
import input_preprocessing.markers_preprocessing as mp

# output csv file
csv_file = '/prodslow/testing/ariel/genomagic_qa/SGU-29/out.csv'

# list of ref fasta files
ref_file_list = ['/prodslow/testing/ariel/genomagic_qa/SGU-29/williams82__ver100.fasta', '/prodslow/testing/ariel/genomagic_qa/SGU-29/lee__ver100.fasta']

# ref names (should be the same length) - should match the ref fasta list, and be the names of the samples in the genomagic
ref_name_list = ['williams82__ver100', 'lee__ver100']

# the markers input fasta in the following format:
# >probe_id_1
# TACGA
# >probe_id_2
# TAGGA
reads_fasta = '/prodslow/testing/ariel/genomagic_qa/SGU-29/nam_6k_design.fasta'

#number of cores to use in the alignment
threads_num = 3

mp.align_to_multiple_genomes_and_write_csv(csv_file, ref_file_list, ref_name_list, reads_fasta, 3)
