#!/usr/bin/python
import sys
simulationFile=sys.argv[1] 
vcfInputFile=sys.argv[2]
vcfOutoutFile=sys.argv[3]
vcfInfoColumnsNum=9
vcfHeaderLinesNum=6

input_sim=open(simulationFile,'r')
input_vcf=open(vcfInputFile,'r')
output_vcf=open(vcfOutoutFile,'w')


line=input_vcf.readline().rstrip()
while line[0:6]!="#CHROM":
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



line1=input_sim.readline().rstrip()
line2=input_vcf.readline().rstrip()
while line1 and line2:
    parts1=line1.split('\t')
    parts2=line2.split('\t')
    start1=int(parts1[1])
    start2=int(parts2[1])
    assert parts1[7][:7]=='HS;END='
    assert parts2[7][:7]=='HS;END='
    end1=int(parts1[7][7:])
    end2=int(parts2[7][7:])
    curr_start=max(start1,start2)
    curr_end=min(end1,end2)
    curr_len=curr_end-curr_start+1
    output_vcf.write('{}\t{}'.format(parts2[0],curr_start))
    for i in range(2,7):
        output_vcf.write('\t{}'.format(parts2[i]))
    output_vcf.write('\tHS;END={}\tGT:CO1:LN'.format(curr_end))
    for i in range(sim_samples_num):
        current_genotype=parts1[vcfInfoColumnsNum+i]
        curr_id=int(current_genotype[:current_genotype.find('|')])-1
        if curr_id<0:
            current_genotype='0|0:0,0,0'
        else:
            current_genotype=parts2[vcfInfoColumnsNum+curr_id]
        output_vcf.write('\t{}:{}'.format(current_genotype[:current_genotype.rfind(':')],curr_len))
    output_vcf.write('\n')
    if end1 == end2:
        line1=input_sim.readline().rstrip()
        line2=input_vcf.readline().rstrip()
    elif end1<end2:
        line1=input_sim.readline().rstrip()
    else:
        line2=input_vcf.readline().rstrip()
        
    
#sim_pos=[-1,-1]
#line=input_vcf.readline().rstrip()
#while line:        
#    parts=line.split('\t')
#    output_vcf.write(parts[0])
#    for i in range(1,vcfInfoColumnsNum):
#        output_vcf.write("\t{}".format(parts[i]))
#    vcf_pos=int(parts[1])
#    while(vcf_pos>sim_pos[1]):        
#        sim_line=input_sim.readline().rstrip().split('\t')
#        assert(len(sim_line)>7),"length of line is {} : {}".format(len(sim_line),sim_line[0])
#        assert(len(sim_line[7])>7)
#        sim_pos[0]=int(sim_line[1])
#        sim_pos[1]=int(sim_line[7][7:])
#    assert len(sim_line)==vcfInfoColumnsNum+sim_samples_num
#    for i in range(sim_samples_num):
#        current_genotype=sim_line[vcfInfoColumnsNum+i]
#        prev_id=int(current_genotype[:current_genotype.find('|')])
#        new_id=prev_id-1#sim_to_vcf_index[prev_id-1]
#        #print "g={}\t prev_id=={} new_id=={}".format(current_genotype,prev_id,new_id)
#        output_vcf.write("\t{}".format(parts[vcfInfoColumnsNum+new_id]))
#    assert sim_pos[0]<=vcf_pos,"sim_pos=[{},{}], vcf_pos={}".format(sim_pos[0],sim_pos[1],vcf_pos)
#    assert sim_pos[1]>=vcf_pos
#    output_vcf.write('\n')
#    line=input_vcf.readline().rstrip()
input_sim.close()
input_vcf.close()
output_vcf.close()

