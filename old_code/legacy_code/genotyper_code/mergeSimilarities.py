#!/usr/bin/python
import sys
import BreedingSimulationUtils
max_chr_size=320000000
similarities_file_name=sys.argv[1]
vcf_file_for_sample_names=sys.argv[2]
simulations_num=int(sys.argv[3]) 
output_filename=sys.argv[4]
f=open(similarities_file_name,'r')
line=f.readline().rstrip()
lines_list=[]
pos_set=set()
pos_set.add(1)
pos_set.add(max_chr_size)
line=f.readline().rstrip()
while line:
    parts=line.split('\t')
    assert len(parts)==4
    for i in range(2,4):
        parts[i]=int(parts[i])
        pos_set.add(parts[i])
    lines_list.append(parts)
    # a = sorted(a, key=lambda x: x.modified, reverse=True)
    line=f.readline().rstrip()
   #line.strsplit
f.close()
lines_list = sorted(lines_list, key=lambda x: x[3]-x[2], reverse=True)
pos_set=sorted(pos_set)
positions_index={}
positions_num=len(pos_set)
for i in range(positions_num):
    positions_index[pos_set[i]]=i
sample_names=GenotyperUtils.extractSampleNamesFromVCF(vcf_file_for_sample_names,9)
initial_group_size=len(sample_names)
sample_index={}
data_matrix = [[0 for x in range(positions_num-1)] for y in range(simulations_num+ initial_group_size)]

for i in range(initial_group_size):
    sample_index[sample_names[i]]=i+1
    data_matrix[i]=[i+1 for x in range(positions_num-1)]

for i in range(simulations_num):
    sample_names.append("sample_{}".format(i+initial_group_size+1))

similarity_lines_num=len(lines_list)
for i in range(similarity_lines_num):
    target_val=sample_index[lines_list[i][0]]
    assert lines_list[i][1][:7]=='sample_'
    curr_sample_index= int(lines_list[i][1][7:])-1
    assert curr_sample_index>=initial_group_size
    assert curr_sample_index<simulations_num+ initial_group_size
    sim_start=positions_index[lines_list[i][2]]
    sim_end=positions_index[lines_list[i][3]]
    for j in range(sim_start,sim_end):
        if data_matrix[curr_sample_index][j]==0:
            data_matrix[curr_sample_index][j]=target_val

BreedingSimulationUtils.writeSimulationsToVcf(output_filename,data_matrix,pos_set,sample_names)


    
    
    
    




