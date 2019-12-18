import os
import sys
working_dir_path = os.getcwd()
genomagic_qa_repo_path = '/'.join(working_dir_path.split('/')[:-1])
sys.path.append(genomagic_qa_repo_path)
import vcf_utils.vcf_iterator as vi

vcf_file = '/prodslow/testing/ariel/genomagic_qa/PSG-16/temp1'
v = vi.VcfIterator(vcf_file)

