#!/usr/bin/env python
import numpy

def selectRandom(data,pos,k):
    n=len(data)
    x=numpy.ceil(n*numpy.random.sample(k)).astype(int)
    return x

def myMin(v):
    n=len(v)
    min_val=v[0]+1
    min_ind=0
    for i in range(n):
        if v[i]<min_val:
            min_val=v[i]
            min_ind=i
    return [min_val,min_ind]


def compare_single_marker(sample,parent1,parent2):
    """"
    if sample==parent1== opposite of parent2 - returns 1
    if sample==parent2== opposite of parent1 - returns 2
    else returns 0
    """
    if sample == 0 or parent1 == 0 or parent2 == 0:
        return 0
    elif parent1 == parent2:
        return 0
    elif sample == parent1:
        return 1
    else:
        assert sample == parent2
        return 2
    
    
def compare_markers(pos, sample_data, parent1_data, parent2_data):
    markers_num = len(pos)
    assert(len(sample_data) == markers_num)
    assert(len(parent1_data) == markers_num)
    assert(len(parent2_data) == markers_num)
    previous_state = 0
    similarities = []
    for i in range(markers_num):
        res = compare_single_marker(sample_data[i],parent1_data[i],parent2_data[i])
        if res > 0:
            if res == previous_state:
                similarities[-1][2] = pos[i]
            else:
                similarities.append([res, pos[i], pos[i]])
            previous_state = res
    return similarities
