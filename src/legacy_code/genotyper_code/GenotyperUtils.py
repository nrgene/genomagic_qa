#!/usr/bin/python
import numpy

def readIntMatrixFromText(input_file):
    matrix=[]
    f = open ( input_file , 'r')
    line=f.readline()
    while line:
        matrix.append(map(int, line.rstrip().split('\t')))     
        line=f.readline()   
    f.close()
    return matrix

def transposeMatrix(matrix):
    rows_num=len(matrix)
    cols_num=len(matrix[0])
    new_matrix = [[0 for x in range(rows_num)] for y in range(cols_num)]
    for i in range(rows_num):
        for j in range(cols_num):
            new_matrix[j][i]=matrix[i][j]
    return new_matrix

def removeDuplicatedLines(matrix,pos):
    rows_num=len(matrix)
    assert(len(pos)==len(matrix)+1)   
    new_matrix=[]
    new_pos=[]
    new_matrix.append(matrix[0][:])
    new_pos.append(pos[0])    
    for i in range(1,rows_num):
        if matrix[i]!=matrix[i-1]:
            new_matrix.append(matrix[i][:])
            new_pos.append(pos[i])
    new_pos.append(pos[-1])        
    return [new_matrix,new_pos]

def mergeSimilarityVCF(vcf_in, vcf_out):
    input_file=open(vcf_in,'r')
    output_file=open(vcf_out,'w')
    line=input_file.readline().rstrip()
    while(line[0]=='#'):
        output_file.write('{}\n'.format(line))
        line=input_file.readline().rstrip()
    #prev_parts=line.split('\t')
    prev_parts=removeLengthFromText(line.split('\t'),9)
    cols_num=len(prev_parts)
    line=input_file.readline().rstrip()
    while line:
        #print prev_parts[9:12]
        parts=removeLengthFromText(line.split('\t'),9)
        #parts=line.split('\t')
        if parts[9:]==prev_parts[9:]:
            #print "im here"
            assert parts[2:7]==prev_parts[2:7]
            parts[1]=prev_parts[1]
        else:
            assert prev_parts[7][:7]=='HS;END='
            bin_length=int(prev_parts[7][7:])-int(prev_parts[1])+1
            output_file.write('{}'.format(prev_parts[0]))
            for i in range(1,9):
                output_file.write('\t{}'.format(prev_parts[i]))
            for i in range(9,cols_num):
                output_file.write('\t{}:{}'.format(prev_parts[i],bin_length))

            output_file.write('\n')
        prev_parts=parts
        line=input_file.readline().rstrip()    
    input_file.close()
    assert prev_parts[7][:7]=='HS;END='
    bin_length=int(prev_parts[7][7:])-int(prev_parts[1])+1
    output_file.write('{}'.format(prev_parts[0]))
    for i in range(1,9):
        output_file.write('\t{}'.format(prev_parts[i]))
    for i in range(9,cols_num):
        output_file.write('\t{}:{}'.format(prev_parts[i],bin_length))

    output_file.write('\n')
    output_file.close()

def removeLengthFromText(arr,start_pos):
    l=len(arr)
    for i in range(start_pos,l):
        a=arr[i]
        arr[i]=a[:a.rfind(':')]
        #arr[i]=a[:a.index(':')]
    return arr

def removeDuplicatedLinesTransposed(matrix,pos):
    tr_matrix=transposeMatrix(matrix)
    [matrix_tr_short, pos_short]=removeDuplicatedLines(tr_matrix,pos)
    matrix_short = transposeMatrix(matrix_tr_short)
    return [matrix_short,pos_short]

def extractSampleNamesFromVCF(filename,infoColumnsNum):
    f = open ( filename , 'r')
    line=f.readline()
    while(line[0]=='#'):
        parts=line.rstrip().split('\t')
        line=f.readline()
    f.close()
    return parts[infoColumnsNum:]

def findSnp(seq1, seq2):
    n1=len(seq1)
    n2=len(seq2)
    i=0
    while seq1[i]==seq2[i]:
        assert i<n1
        assert i<n2
        i+=1
    if (n1==n2):
        assert seq1[i+1:]==seq2[i+1:]
        return [ seq1[i],seq2[i]]
    elif n1>n2:
        d=n1-n2
        assert seq1[i+d:]==seq2[i:]
        return [seq1[i:i+d],'-']
    else:
        d=n2-n1
        assert seq2[i+d:]==seq1[i:]
        return ['-',seq2[i:i+d]]


  
def snpToGenomagicInput(snp_file,tsv_out):
    input_file=open(snp_file,'r')
    output_file=open(tsv_out,'w')
    output_file.write('sample_id\tmarker_id\tallele_seq\tstatus\n')
    line=input_file.readline()
    while line[0]=='#':
        parts=line.rstrip().split('\t')
        line=input_file.readline()
    headers=parts[9:]    
    samples_num=len(headers)
    while line:
        parts=line.rstrip().split('\t')
        seqs=parts[4].split(',')
        [alt1,alt2]=findSnp(seqs[0][1:-1],seqs[1][1:-1])
        if alt1=='-':
            assert len(alt2)>1
            alt1='*'*len(alt2)
        if alt2=='-':
            assert len(alt1)>1
            alt2='*'*len(alt1)
        for i in range(samples_num):            
            assert parts[9+i][0]==parts[9+i][2]
            g=parts[9+i][0]            
            if g=='1':
                output_file.write('{}\t{}\t{}\t1\n'.format(headers[i],parts[2],alt1))
                output_file.write('{}\t{}\t{}\t0\n'.format(headers[i],parts[2],alt2))                
            elif g=='2':
                output_file.write('{}\t{}\t{}\t1\n'.format(headers[i],parts[2],alt2))
                output_file.write('{}\t{}\t{}\t0\n'.format(headers[i],parts[2],alt1))        
            else:
                assert g=='3'        
        line=input_file.readline()
    input_file.close()
    output_file.close()

def estimateRecombinationFromNamSimilarity(nam_file):    
    headers_columns_num=6
    vcf_info_columns_num=11
    resolution=10000000
    chr_size=320000000
    input_file=open(nam_file,'r')
    for i in range(headers_columns_num):
        line=input_file.readline()
    parts=line.rstrip().split('\t')
    samples_num=len(parts)-vcf_info_columns_num
    prev_vec=['0']*samples_num
    line=input_file.readline()
    parts=line.rstrip().split('\t')
    for i in range(samples_num):
        val=parts[vcf_info_columns_num+i][0]
        if val!='3':
            prev_vec[i]=val        
    line=input_file.readline()
    recombination_positions=[];
    while line:
        parts=line.rstrip().split('\t')
        current_pos=int(parts[1])
        for i in range(samples_num):
            val=parts[vcf_info_columns_num+i][0]
            if val!='3':
                if prev_vec[i]=='0':
                    prev_vec[i]=val
                elif prev_vec[i]!=val:
                    assert prev_vec[i]!='0'
                    assert val=='1' or val=='2'                    
                    recombination_positions.append(current_pos)
                    prev_vec[i]=val                
        line=input_file.readline()  
    input_file.close()    
    [h,b]=numpy.histogram(recombination_positions,bins=range(0,chr_size+1,resolution))    
    h = map(lambda x: x /(samples_num*4.0), h)    
    cs=numpy.cumsum(h)
    return [cs,b]
    





