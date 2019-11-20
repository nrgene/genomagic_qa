#!/usr/bin/python

def filterExpression(positions,v1,v2,is_informative):
    n=len(positions)
    assert len(v1)==n
    assert len(v2)==n
    assert len(is_informative)==n    
    new_pos=[]
    similar=[]
    for i in range(n):
        if v1[i]!=4 and v2[i]!=4:# if has_entry_in(reference) AND has_entry_in(compared_to)
            if (is_informative[i][0] and (v1[i]==1 or v2[i]==1)) or (is_informative[i][1] and (v1[i]==2 or v2[i]==2)): # is_informative AND (exists_in(reference) OR exists_in(compared_to))
                new_pos.append(positions[i])
                similar.append(v1[i]==v2[i])
    return [new_pos,similar]

def findSequences(vec,windowSize,maxErrorsNum):
    start_list=[]
    end_list=[]
    n = len(vec)
    slidingErrorCount=0   
    for i in range(windowSize):        
        if not vec[i]:
            slidingErrorCount+=1
    if slidingErrorCount<=maxErrorsNum:
        ind=i-windowSize+1
        while not vec[ind]:
            ind+=1
        start_list.append(ind)
        #print "start: {}".format(ind)
    for i in range(windowSize,n):
        slidingErrorCountPrev=slidingErrorCount        
        if vec[i]:            
            if not vec[i-windowSize]:
                slidingErrorCount-=1
        else:
            assert not vec[i]
            if vec[i-windowSize]:
                slidingErrorCount+=1
        assert slidingErrorCount>=0,"i={}\tslidingErrorCount={}".format(i,slidingErrorCount)
        
        if slidingErrorCount<=maxErrorsNum and slidingErrorCountPrev>maxErrorsNum:
            ind=i-windowSize+1
            while not vec[ind]:
                ind+=1
            start_list.append(ind)
            #print "start: {}".format(ind)
        elif slidingErrorCount>maxErrorsNum and slidingErrorCountPrev<=maxErrorsNum:
            ind=i-1
            while not vec[ind]:
                ind-=1
            end_list.append(ind)
            #print "end: {}".format(ind)       
            
    if slidingErrorCount<=maxErrorsNum:
        ind=n-1
        while not vec[ind]:
            ind-=1
        end_list.append(ind)
    return [start_list,end_list]
        #print "end: {}".format(ind)
    
    

def snpSimilarity(positions,v1,v2,is_informative,windowSize,maxErrorsNum):
    [new_pos,similar]=filterExpression(positions,v1,v2,is_informative)
    [start_list,end_list]=findSequences(similar,windowSize,maxErrorsNum)
    m=len(start_list)
    assert len(end_list)==m
    similarities=[]    
    for i in range(m):
        start=new_pos[start_list[i]]
        end=new_pos[end_list[i]]
        assert start<=end,"start={} >= end={}".format(start,end)
        similarities.append((start,end))
    return similarities
        
        
def readSnpFromFile(filename):
    input_sim=open(filename,'r')
    sim_line=input_sim.readline().rstrip()
    while sim_line[0:6]!="#CHROM":
        sim_line=input_sim.readline().rstrip()
    parts=sim_line.split('\t')
    sim_samples=parts[9:]
    samples_num= len(sim_samples)
    dataMatrix = [[] for y in range(samples_num)] 
    is_informative=[]
    positions=[]
    sim_line=input_sim.readline().rstrip()
    while sim_line:
        parts=sim_line.split('\t')
        sim_line=input_sim.readline().rstrip()
        #print sim_line
        for i in range(samples_num):
            
            dataMatrix[i].append(int(parts[9+i][0]))
            info=parts[7].split(';')
            info_tuple=(int(info[1][-1]),int(info[2][-1]))
            assert info_tuple[0]==1 or info_tuple[0]==0
            assert info_tuple[1]==1 or info_tuple[1]==0
        is_informative.append(info_tuple)
        positions.append(int(parts[1]))
    input_sim.close()
    return [dataMatrix,positions,is_informative,sim_samples]


    

    
import sys
inputFile=sys.argv[1]
outputFile=sys.argv[2]
group1_start=int(sys.argv[3])
group1_end=int(sys.argv[4])
group2_start=int(sys.argv[5])
group2_end=int(sys.argv[6])


[dataMatrix,positions,is_informative,sim_samples]=readSnpFromFile(inputFile)
input_sim=open(outputFile,'w')
group1=range(group1_start,group1_end)
group2=range(group2_start,group2_end)
n1=len(group1)
n2=len(group2)
for i in group1:
    for j in group2:
        if i!=j:
            similarities=snpSimilarity(positions,dataMatrix[i],dataMatrix[j],is_informative,20,0)
            n=len(similarities)
            for k in range(n):
                input_sim.write("{}\t{}\t{}\t{}\n".format(sim_samples[i],sim_samples[j],similarities[k][0],similarities[k][1]))
input_sim.close()      

