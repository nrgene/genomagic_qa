#!/usr/bin/python
#import sys,re,glob,os
#import re,SNPUtils

def findPermutationFromAtoB(samplesA, samplesB):
    lengthOfA=len(samplesA)
    indexFromAtoB=[-1]*lengthOfA
    for i in range(lengthOfA):
        assert samplesA[i] in samplesB
        indexFromAtoB[i]=samplesB.index(samplesA[i])
        #print "{} was in {} now in {}".format(samplesA[i],i,indexFromAtoB[i])
    assert not -1 in indexFromAtoB
    return indexFromAtoB
    
    

import sys
simulationFile=sys.argv[1] # vcf simulations file 
vcfInputFile=sys.argv[2]#'snp_markers.vcf'
vcfOutoutFile=sys.argv[3]#'temp'
#number_of_initial_samples=int(sys.argv[4])
vcfInfoColumnsNum=9
vcfHeaderLinesNum=6

#vcf_info_rows_num=int(sys.argv[4]) # should be 6 in snp and 7 in hap

input_sim=open(simulationFile,'r')
input_vcf=open(vcfInputFile,'r')
output_vcf=open(vcfOutoutFile,'w')


line=input_vcf.readline().rstrip()
while line[0:vcfHeaderLinesNum]!="#CHROM":
    output_vcf.write('{}\n'.format(line))
    parts=line.split('\t')
    line=input_vcf.readline().rstrip()
parts=line.split('\t')
vcf_samples=parts[vcfInfoColumnsNum:]
sourceSamplesNum=len(vcf_samples)

sim_line=input_sim.readline().rstrip()
while sim_line[0:vcfHeaderLinesNum]!="#CHROM":
    parts=sim_line.split('\t')
    sim_line=input_sim.readline().rstrip()
output_vcf.write('{}\n'.format(sim_line))
parts=sim_line.split('\t')
sim_samples=parts[vcfInfoColumnsNum:]
sim_samples_num=len(sim_samples)

# this is needed when we want to compare sample names
#sim_to_vcf_index=findPermutationFromAtoB(sim_samples[:sourceSamplesNum],vcf_samples)

sim_pos=[-1,-1]
line=input_vcf.readline().rstrip()
while line:        
    parts=line.split('\t')
    output_vcf.write(parts[0])
    for i in range(1,vcfInfoColumnsNum):
        output_vcf.write("\t{}".format(parts[i]))
    vcf_pos=int(parts[1])
    while(vcf_pos>sim_pos[1]):        
        sim_line=input_sim.readline().rstrip().split('\t')
        assert(len(sim_line)>7),"length of line is {} : {}".format(len(sim_line),sim_line[0])
        assert(len(sim_line[7])>7)
        sim_pos[0]=int(sim_line[1])
        sim_pos[1]=int(sim_line[7][7:])
        haps_id=sim_line[4][4:-9].split('>,<HAP')        
    assert len(sim_line)==vcfInfoColumnsNum+sim_samples_num
    for i in range(sim_samples_num):
        current_genotype=sim_line[vcfInfoColumnsNum+i]
        prev_id=int(current_genotype[:current_genotype.find('|')])        
        if prev_id==0 or prev_id>sourceSamplesNum:
            output_vcf.write("\t3|3:0,0,0")
        else:
            new_id=int(haps_id[prev_id-1])-1
            #new_id=prev_id-1#sim_to_vcf_index[prev_id-1]            
            output_vcf.write("\t{}".format(parts[vcfInfoColumnsNum+new_id]))
    assert sim_pos[0]<=vcf_pos,"sim_pos=[{},{}], vcf_pos={}".format(sim_pos[0],sim_pos[1],vcf_pos)
    assert sim_pos[1]>=vcf_pos
    output_vcf.write('\n')
    line=input_vcf.readline().rstrip()
input_sim.close()
input_vcf.close()
output_vcf.close()

