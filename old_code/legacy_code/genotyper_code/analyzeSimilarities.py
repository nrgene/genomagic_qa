#!/usr/bin/python
import VcfUtils
import sys

def computeSimilaritiesBetweenTwoSamples(snp_marker_values,ind1,ind2,window_length,p):
    markers_num=len(snp_marker_values)
    vec=[]
    indices=[]
    i=0
    while i<markers_num:
        #print '{}: ({},{})'.format(i,snp_marker_values[i][ind1],snp_marker_values[i][ind2])
        if snp_marker_values[i][ind1]!=0 and snp_marker_values[i][ind2]!=0:
            vec.append(int(snp_marker_values[i][ind1]==snp_marker_values[i][ind2]))
            indices.append(i)
        i+=1
    #print 'initial length={} final length={}'.format(markers_num,len(vec))
    #print vec
    markers_num=len(vec)
    curr_sum=sum(vec[0:window_length])
    curr_state = curr_sum>=window_length*p
    similarities=[]
    if curr_state:        
        similarities.append([indices[0],-1])
    for i in range(window_length,markers_num):        
        curr_sum+=vec[i]
        curr_sum-=vec[i-window_length]
        prev_state = curr_state
        curr_state = curr_sum>=window_length*p
        if (not prev_state) and curr_state:
            similarities.append([indices[i-window_length+1],-1])
        elif (not curr_state) and prev_state:
            similarities[-1][1]=indices[i-1]
    if curr_state:
        similarities[-1][1]=indices[markers_num-1]    
    return similarities
        
def computeSimilaritiesInSingleChromosome(curr_chr,snp_marker_values,pos,window_length,p):
    
    samples_num=len(snp_marker_values[0])
    comparisons_num=samples_num*(samples_num-1)/2
    total_similarities=[]
    for i in range(samples_num-1):
        for j in range(i+1,samples_num):
            #print 'comp: {},{}'.format(i,j)
            curr_similarities=computeSimilaritiesBetweenTwoSamples(snp_marker_values,i,j,window_length,p)
            for curr_sim in curr_similarities:                
                total_similarities.append([i,j]+[curr_chr,pos[curr_sim[0]][1],pos[curr_sim[1]][2]])
    return total_similarities
            
def getIndicesOfChromosomes(pos):
    n=len(pos)
    prev_chr=pos[0][0]
    prev_pos=pos[0][1]
    X=[]
    X.append([prev_chr,0,-1])
    for i in range(n):
        if pos[i][0]!=prev_chr:
            assert pos[i][0]>=prev_chr            
            X[-1][2]=i-1
            X.append([pos[i][0],i,-1])
        else:
            assert pos[i][1]>=prev_pos
        prev_chr=pos[i][0]
        prev_pos=pos[i][1]
    X[-1][2]=n-1
    for x in X:
        vec=pos[x[1]:(x[2]+1)]
        for v in vec:
            assert v[0]==x[0]
    return X
            
    
def computeSimilaritiesPerEachChromosome(data,pos,sample_names,window_length,p):
    similarity_list=[]
    X=getIndicesOfChromosomes(pos)
    for x in X:
        curr_chr=x[0]
        curr_sim_list=computeSimilaritiesInSingleChromosome(curr_chr,data[x[1]:x[2]+1],pos[x[1]:x[2]+1],window_length,p)
        for sim_record in curr_sim_list:
            similarity_list.append([sample_names[sim_record[0]],sample_names[sim_record[1]]]+sim_record[2:])
    return similarity_list

    
snp_file_name=sys.argv[1]
out_file_name=sys.argv[2]
window_length=int(sys.argv[3])
p=float(sys.argv[4])
[data,pos,informative,sample_names] = VcfUtils.readSnpVcfFile(snp_file_name)
similarity_list=computeSimilaritiesPerEachChromosome(data,pos,sample_names,window_length,p)
f = open(out_file_name, 'w')
for s in similarity_list:
    f.write('{}\t{}\t{}\t{}\t{}\n'.format(s[0],s[1],s[2],s[3],s[4]))
f.close()


    

