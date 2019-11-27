import numpy as np
import io
import os
import pandas as pd


def read_snp(str):
    if str[:3] == '1|1':
        return 1
    elif str[:3] == '2|2':
        return 2
    elif str[:3] == '1|2' or str[:3]=='2|1':
        return 3
    elif str[:3] == '3|3' or str[:3]=='4|4':
        return 0
    else:
        raise NameError(str[:3] + 'not recognized')


def read_vcf_data(vcf_line):
    nsamples = len(vcf_line)
    var_data = np.zeros([1, nsamples], dtype=np.int8)
    for i in range(nsamples):
        var_data[0][i] = read_snp(vcf_line[i])
    return var_data


def read_snp_markers_vcf_file(markers_file, max_lines_number):
    with open(markers_file) as fp:
        prev_line = ''
        line = fp.readline().rstrip()
        while line[0] == '#':
            prev_line = line
            line = fp.readline().rstrip()
        headers_parts = prev_line.split('\t')
        samples = headers_parts[9:]
        nsamples = len(samples)
        data = np.zeros([max_lines_number, nsamples], dtype=np.int8)
        pos = np.zeros([max_lines_number, 2], dtype=np.int)
        count = 0
        while line:
            parts = line.split('\t')
            data[count] = read_vcf_data(parts[9:])
            pos[count] = [int(parts[0]), int(parts[1])]
            line = fp.readline().rstrip()
            count += 1
            if (count % 10000) == 0:
                print("read {:d}K lines".format(int(count/1000)))
            assert count < max_lines_number
        return [pos[:count], data[:count], samples]


def read_file_to_dataframe(markers_file, max_lines_number):
    [pos, data, samples] = read_snp_markers_vcf_file(markers_file, max_lines_number)
    nsamples = len(samples)
    vcf_dict = {'chr': pos[:,0], 'pos' : pos[:,1]}
    for i in range(nsamples):
        vcf_dict[samples[i]] = data[:,i]
        df = pd.DataFrame(vcf_dict)
    return df

