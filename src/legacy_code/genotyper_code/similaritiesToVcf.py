#!/usr/bin/python
import sys
simulationFile=sys.argv[1]
similarity_file=open(simulationFile,'r')
lines_list=[]
position_set=set()
line=similarity_file.readline()
while line:
    parts=line.rstrip().split('\t')
    start=int(parts[2])
    end=int(parts[3])
    position_set.add(start)
    position_set.add(end)    
    lines_list.append((parts[0],parts[1],start,end,end-start))
    line=similarity_file.readline()
   


similarity_file.close()


