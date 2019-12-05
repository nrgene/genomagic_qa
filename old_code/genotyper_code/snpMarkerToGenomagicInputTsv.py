#!/usr/bin/python
import sys
import GenotyperUtils
snp_file=sys.argv[1]
tsv_out=sys.argv[2]
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
    [alt1,alt2]=GenotyperUtils.findSnp(seqs[0][1:-1],seqs[1][1:-1])
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
    #print '{}\t{}'.format(alt1,alt2)
    line=input_file.readline()
input_file.close()
output_file.close()

