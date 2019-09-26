from sets import Set
number_of_info_columns=9
	
def checkMarkerInBin(sample_has_marker,sample_similarity_id):   
    samples_num=len(sample_has_marker)
    assert len(sample_similarity_id)==samples_num
    exists=Set()
    not_exists=Set()
    for k in range(samples_num):
        assert sample_similarity_id[k]>0
        if sample_has_marker[k]:
            exists.add(sample_similarity_id[k])
        else:
            not_exists.add(sample_similarity_id[k])
    exists=list(exists)
    not_exists=list(not_exists)
    #print 'len exist:{} len not exist{}'.format(len(exists),len(not_exists))
    assert len(exists)>0
    if len(exists)>1:
        return -1
    elif exists[0] in not_exists:
        return -2
    else:
        return exists[0]
		
		
def parseHapMarkerLine(line):
    assert line
    parts=line.split('\t')
    chr_id=int(parts[0])
    pos=int(parts[1])
    k=len(parts)-9
    has_marker=[False]*k
    for i in range(k):
        if parts[9+i][0]=='2':
            has_marker[i]=True
        else:
            assert parts[9+i][0]=='1'
    return [chr_id,pos,has_marker]

	
def skipMarkersThatAreBeforeSimilarityBin(similarity_pos,marker_line,markers_file):
    assert marker_line
    [marker_chr,marker_pos,sample_has_marker]=parseHapMarkerLine(marker_line)    
    while marker_line and marker_chr<similarity_pos[0]:        
        marker_line=markers_file.readline().rstrip()
        if marker_line:
            [marker_chr,marker_pos,sample_has_marker]=parseHapMarkerLine(marker_line)
    while marker_line and marker_chr==similarity_pos[0] and marker_pos<similarity_pos[1]:        
        marker_line=markers_file.readline().rstrip()
        if marker_line:
            [marker_chr,marker_pos,sample_has_marker]=parseHapMarkerLine(marker_line)
    if not marker_line:
        print "skipMarkersThatAreBeforeSimilarityBin reached to end of file, but couldnt find markers that are not before {}".format(similarity_pos)    
    return marker_line
	
def readVcfHeader(f,fout):
    number_of_info_columns=9
    line=f.readline().rstrip()
    while line[0]=="#":
        fout.write('{}\n'.format(line))
        parts=line.split('\t')
        line=f.readline().rstrip()
    info_fields=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
    assert parts[:number_of_info_columns]==info_fields,parts[:number_of_info_columns]
    return [line,parts]
	

def handleCurrentSimilarityBin(similarity_data,similarity_pos,marker_line,markers_file,out_file,log_file):    
    hap_set=Set(similarity_data)
    total_haps_in_bin=len(hap_set)
    if marker_line:  
        marker_line=skipMarkersThatAreBeforeSimilarityBin(similarity_pos,marker_line,markers_file)    
    if marker_line:        
        [marker_chr,marker_pos,sample_has_marker]=parseHapMarkerLine(marker_line)
        assert marker_chr>=similarity_pos[0]
        assert marker_pos>=similarity_pos[1]
        while marker_line and marker_chr==similarity_pos[0] and marker_pos<similarity_pos[2]:         
            assert marker_pos>=similarity_pos[1]
            marker_info=checkMarkerInBin(sample_has_marker,similarity_data)            
            if marker_info in hap_set:               
                hap_set.remove(marker_info)            
                log_file.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(similarity_pos[0],similarity_pos[1],similarity_pos[2],marker_info,marker_chr,marker_pos))
                out_file.write('{}\n'.format(marker_line))
            marker_line=markers_file.readline().rstrip()
            if marker_line:
                [marker_chr,marker_pos,sample_has_marker]=parseHapMarkerLine(marker_line)
    #else:
        #print "no markers found for similarity {}:{}-{}".format(similarity_pos[0],similarity_pos[1],similarity_pos[2])
    for k in hap_set:
        log_file.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(similarity_pos[0],similarity_pos[1],similarity_pos[2],k,-1,-1))                          
    unmatched_haps_in_bin=len(hap_set)
    return [unmatched_haps_in_bin,total_haps_in_bin,marker_line]

def iterateOverSimilarities(similarity_data,similarity_pos,sample_names,markers_file,output,logfile):
    samples_num=len(sample_names)
    similarity_bins_num=len(similarity_data)
    assert len(similarity_pos)==similarity_bins_num
    assert len(similarity_data[0])==samples_num
    [marker_line,parts]=readVcfHeader(markers_file,output)
    assert parts[9:]==sample_names,'{}\n{}'.format(parts[9:],sample_names)    
    total_unmatched=0
    total_haps=0
    for i in range(similarity_bins_num):        
        [unmatched_haps_in_bin,total_haps_in_bin,marker_line]=handleCurrentSimilarityBin(similarity_data[i],similarity_pos[i],marker_line,markers_file,output,logfile)        
        total_unmatched+=unmatched_haps_in_bin    
        total_haps+=total_haps_in_bin
    return [total_unmatched,total_haps]  



def readVcfHeader1(f):    
    line=f.readline().rstrip()
    while line[0]=="#":
        parts=line.split('\t')
        line=f.readline().rstrip()
    info_fields=["#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT"]
    assert parts[:number_of_info_columns]==info_fields,parts[:number_of_info_columns]
    return [line,parts]

def readHaplotypeSimilarityVcfFile(f):    
    [line,parts]=readVcfHeader1(f)    
    sample_names=parts[number_of_info_columns:]
    samples_num=len(sample_names)
    data=[]
    pos=[]    
    while line:
        parts=line.split('\t')
        b=parts[7].split('=')
        pos.append([int(parts[0]),int(parts[1]),int(b[1])])
        arr=[0]*samples_num
        for i in range(samples_num):
            b=parts[number_of_info_columns+i].split('|')
            arr[i]=int(b[0])
        data.append(arr)   
        line=f.readline().rstrip()  
    f.close()
    return [data,pos,sample_names]
