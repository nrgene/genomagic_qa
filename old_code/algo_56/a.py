#!/usr/bin/env python
import markersSelection
import sys
print "hello"
sys.path.append('C:\\Users\\ariel\\workspace\\arielUtils')
import VcfUtils

filename = 'wgs_snp_markers.vcf'
[data, pos, informative, sample_names] = VcfUtils.readSnpVcfFile(filename)

p = [(row[1] + row[2]) / 2 for row in pos[0:1500]]
sample_data = [row[0] for row in data[0:1500]]
parent1_data = [row[1] for row in data[0:1500]]
parent2_data = [row[2] for row in data[0:1500]]
similarities = markersSelection.compareMarkers(p, sample_data, parent1_data, parent2_data)
print similarities
# K=500
# all_markers_num=len(data)
# markers_filtered=[]
# filtered_pos=[]
# parents=[4,17]


# I=markersSelection.selectRandom(data,pos,K)
# for i in range(all_markers_num):
#    if np.mean(data[i]==1)<0.7 and np.mean(data[i]==-1)<0.7 and sum(data[i]==0)==0:
#        markers_filtered.append(data[i])
#        filtered_pos.append(np.mean(pos[i][1:]))
# X=abs(np.cov(markers_filtered))
# n=len(X)
# i=int(np.ceil(len(X)*np.random.random_sample()))
# hap_indices=[i]
# while len(hap_indices)<K:
#    p=[0]*n
#    for i in range(n):







