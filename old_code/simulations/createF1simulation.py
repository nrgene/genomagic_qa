#!/usr/bin/python
import sys,os,random
sys.path.append('/prodslow/testing/ariel/genotyper/code')
import BreedingSimulationUtils,GenotyperUtils

def createCrossingSimulations(nam_file,sample_names):    
    chr_size=320000001
    #crossings_num=400
    samples_per_crossing=50
    initial_group_size=len(sample_names)
    #initial_group_size=400
    total_samples_num=initial_group_size*(samples_per_crossing+1)
    new_sample_names=[]
    recombinations=BreedingSimulationUtils.returnRandomRecombinations(nam_file,initial_group_size*samples_per_crossing)
    pos_set=set()
    pos_set.update(recombinations)
    pos_set.remove(0)
    pos_set.update([1,chr_size])
    pos_set=sorted(pos_set)
    positions_num=len(pos_set)-1
    positions_index = {}
    for i in range(positions_num):
        positions_index[pos_set[i]]=i
    data_matrix = [[0 for x in range(positions_num)] for y in range(total_samples_num)] 
    for i in range(initial_group_size):
        data_matrix[i]=[i+1]*positions_num
        #new_sample_names.append(sample_names[i])
    parents_pool=range(initial_group_size)+[0]
    for current_crossing in range(initial_group_size):
        parents=parents_pool[current_crossing:current_crossing+2]
        for current_sample_in_crossing in range(samples_per_crossing):
            new_sample_names.append('{}_{}_{}'.format(sample_names[parents[0]],sample_names[parents[1]],current_sample_in_crossing+1))
            curr_index=samples_per_crossing*current_crossing+current_sample_in_crossing                
            r=recombinations[curr_index]
            p_vec=parents[:]
            random.shuffle(p_vec)
            data_matrix[curr_index+initial_group_size]=data_matrix[p_vec[0]][:]
            if r>0:
                curr_ind=positions_index[r]
                data_matrix[curr_index+initial_group_size][curr_ind:]=data_matrix[p_vec[1]][curr_ind:]
    data_matrix=data_matrix[initial_group_size:]
    assert len(data_matrix)==len(new_sample_names)
    return [data_matrix,pos_set,new_sample_names]

nam_file='/prodslow/testing/ariel/genotyper/ALGO-39/nam_similarity.vcf'
#hs_file='/prodslow/testing/ariel/genotyper/ALGO-39/haplotype_similarity.vcf'
founders_file='/prodslow/testing/ariel/genotyper/ALGO-45/simulations_of_founder_lines.vcf'
out_file='/prodslow/testing/ariel/genotyper/ALGO-45/simulations_of_F1_lines.vcf'
chr_size=320000001
vcf_file_for_sample_names=founders_file
sample_names=GenotyperUtils.extractSampleNamesFromVCF(vcf_file_for_sample_names,9)
initial_group_size=len(sample_names)
[data_matrix,pos_set,sample_names]=createCrossingSimulations(nam_file,sample_names)#BreedingSimulationUtils.createSimulations(nam_file,chr_size)
#data_matrix_len=len(data_matrix)
#for i in range(initial_group_size,data_matrix_len):
#    sample_names.append('sample_{}'.format(i+1))
BreedingSimulationUtils.writeSimulationsToVcf(out_file,data_matrix,pos_set,sample_names)
