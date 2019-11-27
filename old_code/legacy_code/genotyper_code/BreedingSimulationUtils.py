#!/usr/bin/python
import numpy,sys
import GenotyperUtils

def createSimulations(nam_file,chr_size,initial_group_size,generations_num,samples_per_generation):    
    total_samples_num=initial_group_size+generations_num*samples_per_generation
    recombinations=returnRandomRecombinations(nam_file,generations_num*samples_per_generation)
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
    parents_pool=[0,initial_group_size]
    for curr_gen in range(generations_num):
        for current_sample_in_generation in range(samples_per_generation):
            curr_index=samples_per_generation*curr_gen+current_sample_in_generation
            parents=numpy.random.randint(parents_pool[0],parents_pool[1],2)
            r=recombinations[curr_index]
            data_matrix[curr_index+initial_group_size]=data_matrix[parents[0]][:]
            if r>0:
                assert(r in positions_index),'{} is not in positions_index dict'.format(len(positions_index))
                curr_ind=positions_index[r]
                data_matrix[curr_index+initial_group_size][curr_ind:]=data_matrix[parents[1]][curr_ind:]                
        parents_pool=[initial_group_size+samples_per_generation*curr_gen,initial_group_size+samples_per_generation*(curr_gen+1)]                   
    return [data_matrix,pos_set]
   
def returnRandomRecombinations(nam_file,n):
    [cs,b]=GenotyperUtils.estimateRecombinationFromNamSimilarity(nam_file)
    resolution=1
    m=len(cs)
    I=numpy.digitize(numpy.random.uniform(0,1,n),cs)    
    for i in range(n):
        if I[i]==m:
            I[i]=0
        else:
            I[i]=numpy.random.randint(b[I[i]],b[I[i]+1])
            I[i]=(I[i]/resolution)*resolution
    return I

def writeSimulationsToVcf(filename,matrix,pos,samples_name):    
    [simulations_matrix,positions]=GenotyperUtils.removeDuplicatedLinesTransposed(matrix,pos)
    #colorset=GenotyperUtils.readIntMatrixFromText('/cygdrive/c/Users/ariel/workspace/ALGO-39/NRGeneColorsSet.tsv') 
    colorset=GenotyperUtils.readIntMatrixFromText('/home/ariel/clients/genomagicapi/apicore/src/main/resources/NRGeneColorsSet.tsv')
    bins_num=len(positions)-1    
    samples_num=len(simulations_matrix)
    assert len(samples_name)==samples_num,'len(samples_name)[{}] ==samples_num[{}]'.format(len(samples_name),samples_num)
    assert bins_num==len(simulations_matrix[0])
    max_val=max(max(simulations_matrix))
    hap_text=''
    for i in range(max_val):
        hap_text=hap_text+"<HAP{}>,".format(i+1)
    hap_text=hap_text+"<NAHAP>"
    f = open ( filename , 'w')
    f.write('##fileformat=VCFv4.0\n')
    f.write('##FORMAT=<ID=CO1,Number=3,Type=Integer,Description="Haplotypes Similarity Based Coloring (RGB)">\n')
    f.write('##FORMAT=<ID=LN,Number=1,Type=Integer,Description="Overrides the length (in bases). Optional. It can be present in some or all the variants on the line">\n')
    f.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">\n')
    f.write('##INFO=<ID=HS,Number=0,Type=Flag,Description="Haplotype similarity variant type">\n')
    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT')
    for i in range( samples_num):
        f.write('\t{}'.format(samples_name[i]))
    f.write('\n')
    for i in range(bins_num):
        l=positions[i+1]-positions[i]
        f.write("1\t{}\t.\tN\t{}\t1\tPASS\tHS;END={}\tGT:CO1:LN".format(positions[i],hap_text,positions[i+1]-1))
        for j in range(samples_num):
            val=simulations_matrix[j][i]            
            if val==0:
                f.write("\t{}|{}:0,0,0:{}".format(val,val,l))                
            else:
                f.write("\t{}|{}:{},{},{}:{}".format(val,val,colorset[val-1][0],colorset[val-1][1],colorset[val-1][2],l))                
        f.write('\n')
    f.close()

