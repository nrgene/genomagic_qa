#!/usr/bin/python
import sys,os
sys.path.append('/prodslow/testing/ariel/genotyper/code')
import BreedingSimulationUtils,GenotyperUtils


nam_file='/prodslow/testing/ariel/genotyper/ALGO-39/nam_similarity.vcf'
hs_file='/prodslow/testing/ariel/genotyper/ALGO-39/haplotype_similarity.vcf'
out_file='/prodslow/testing/ariel/genotyper/ALGO-45/simulations_of_founder_lines.vcf'
chr_size=320000001
generations_num=10
samples_per_generation=100
vcf_file_for_sample_names=hs_file
sample_names=GenotyperUtils.extractSampleNamesFromVCF(vcf_file_for_sample_names,9)
initial_group_size=len(sample_names)
[data_matrix,pos_set]=BreedingSimulationUtils.createSimulations(nam_file,chr_size,initial_group_size,generations_num,samples_per_generation)
data_matrix=data_matrix[0:34]+data_matrix[668:1034]
data_matrix_len=len(data_matrix)
sample_names=[]
for i in range(0,data_matrix_len):
    sample_names.append('f{}'.format(i+1))
BreedingSimulationUtils.writeSimulationsToVcf(out_file,data_matrix,pos_set,sample_names)
