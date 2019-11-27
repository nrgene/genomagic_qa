#!/usr/bin/python


def calculateNumberOfSimilarities(parts1,parts2):
    matching=0
    different=0
    n=len(parts1)
    assert len(parts2)==n
    for i in range(n-9):
        x1=parts1[9+i].split('|')
        x2=parts2[9+i].split('|')
        if x1[0]==x2[0]:
            matching+=1
        else:
            different+=1
    return [matching,different]

def readNextLines(pos1,pos2,input1,input2,sim1_line,sim2_line):
    if (pos1[0]<pos2[0]) or (pos1[0]==pos2[0] and pos1[2]<pos2[2]):
        sim1_line=input1.readline().rstrip()
    elif (pos2[0]<pos1[0]) or (pos1[0]==pos2[0] and pos2[2]<pos1[2]):
        sim2_line=input2.readline().rstrip()
    else:
        assert pos2[0]==pos1[0] and pos2[2]==pos1[2]
        sim1_line=input1.readline().rstrip()
        sim2_line=input2.readline().rstrip()
    return [sim1_line,sim2_line]
    


def compareHaplotypeSimilarityFiles(input1,input2):    
    sim1_line=input1.readline().rstrip()
    sim2_line=input2.readline().rstrip()
    while(sim1_line[0]=='#'):    
        sim1_line=input1.readline().rstrip()
    while(sim2_line[0]=='#'):    
        sim2_line=input2.readline().rstrip()
    parts1=sim1_line.split('\t')
    parts2=sim2_line.split('\t')
    pos1=[int(parts1[0]),int(parts1[1]),int(parts1[7][7:])]
    pos2=[int(parts2[0]),int(parts2[1]),int(parts2[7][7:])]
    matching=0
    different=0
    while sim1_line and sim2_line:
        parts1=sim1_line.split('\t')
        parts2=sim2_line.split('\t')
        pos1=[int(parts1[0]),int(parts1[1]),int(parts1[7][7:])]
        pos2=[int(parts2[0]),int(parts2[1]),int(parts2[7][7:])]
        assert pos1[0]==pos2[0]
        assert pos1[1]<pos2[2]
        assert pos2[1]<pos1[2]    
        start_pos=max(pos1[1],pos2[1])
        end_pos=min(pos1[2],pos2[2])
        assert start_pos<=end_pos
        segment_length=end_pos-start_pos+1
        [m,d]=calculateNumberOfSimilarities(parts1,parts2)
        matching+=segment_length*m
        different+=segment_length*d
        [sim1_line,sim2_line]=readNextLines(pos1,pos2,input1,input2,sim1_line,sim2_line)    
    total_length=float(matching+different)
    return round(matching/total_length,2)

#atching % is {}'.format(round(matching/total_length,2))

#similarityFile1='similarity_snp.vcf'
#similarityFile2='similarity_wgs.vcf'
#input1=open(similarityFile1,'r')
#input2=open(similarityFile2,'r')
#r=compareHaplotypeSimilarityFiles(input1,input2)
#input1.close()
#input2.close()
#print 'matching % is {}'.format(r)
