#!/usr/bin/python

vcf1=sys.argv[1]
vcf2=sys.argv[2]
#'similarity_first_10.vcf';
#vcf2='simulations_of_F1_lines.vcf';
#ind1=21:520;
#ind2=10:509;
#samples_num=length(ind1);
#assert(length(ind2)==samples_num);
#Z=zeros(3,samples_num);

#fid1=fopen(vcf1,'r');
#fid2=fopen(vcf2,'r');
#for i = 1:6
#    line1=fgetl(fid1);
#    line2=fgetl(fid2);
#end
#parts1=strsplit(line1,'\t');
#parts2=strsplit(line2,'\t');
#assert(min(strcmp(parts1(ind1),parts2(ind2)))==1);

#line1=fgetl(fid1);
#line2=fgetl(fid2);
#while and(ischar(line1),ischar(line2))    
#    parts1=strsplit(line1,'\t');
#    parts2=strsplit(line2,'\t');
#    start1=str2double(parts1{2});
#    start2=str2double(parts2{2});
#    end1=str2double(parts1{8}(8:end));
#    end2=str2double(parts2{8}(8:end));    
#    curr_len=min(end1,end2)-max(start1,start2);
#    assert(curr_len>0)
#    for i = 1:samples_num
#        if regexp(parts1{ind1(i)},':0,0,0:')
#            Z(1,i)=Z(1,i)+curr_len;
#        else
#            allele1=parts1{ind1(i)}(1:regexp(parts1{ind1(i)},'\|')-1);
#            allele2=parts2{ind2(i)}(1:regexp(parts2{ind2(i)},'\|')-1);
#            if strcmp(allele1,allele2)
#                Z(2,i)=Z(2,i)+curr_len;
#            else
#                Z(3,i)=Z(3,i)+curr_len;
#            end
#        end
#    end    
#    if end1 == end2
#        disp(end1)
#        line1=fgetl(fid1);
#        line2=fgetl(fid2);        
#    elseif end1<end2
#        disp(end1)
#        line1=fgetl(fid1);        
#    else
#        line2=fgetl(fid2);        
#    end   
#end
#fclose(fid1);
#fclose(fid2);
