#!/usr/bin/python
import numpy,sys
import GenotyperUtils
import BreedingSimulationUtils


nam_file='../ALGO-39/nam_similarity.vcf'
chr_size=320000000
generations_num=int(sys.argv[1])
samples_per_generation=int(sys.argv[2])
vcf_file_for_sample_names=sys.argv[4]
sample_names=GenotyperUtils.extractSampleNamesFromVCF(vcf_file_for_sample_names,9)
#'../ALGO-39/haplotype_similarity.vcf',9)
initial_group_size=len(sample_names)
out_file=sys.argv[3]
[data_matrix,pos_set]=BreedingSimulationUtils.createSimulations(nam_file,chr_size,initial_group_size,generations_num,samples_per_generation)
for i in range(initial_group_size,initial_group_size+generations_num*samples_per_generation):
    sample_names.append('sample_{}'.format(i+1))
BreedingSimulationUtils.writeSimulationsToVcf(out_file,data_matrix,pos_set,sample_names)
