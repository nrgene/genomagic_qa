import numpy as np
import sys

def genotype_is_imputed(genotype):
    return genotype.split(':')[1].split(',')[0] == '170'


def handle_new_markers(test_line, stats):
    parts = test_line.split('\t')
    samples_num = len(parts) - 9
    for i in range(samples_num):
        if genotype_is_imputed(parts[i + 9]):
            stats[i][1] += 1


def handle_missing_markers(ref_line, stats):
    parts = ref_line.split('\t')
    samples_num = len(parts)-9
    for i in range(samples_num):
        if parts[i + 9][0] == '1' or parts[i + 9][0] == '2':
            stats[i][0] += 1
            stats[i][4] += 1


def compare_snp_marker_in_sample(expected, observed):
    is_imputed = genotype_is_imputed(observed)
    is_missing = observed[0] == '3' or observed[0] == '4'
    has_ref = expected[0] == '1' or expected[0] == '2'
    if is_imputed:
        if not has_ref:
            return 1  # new value
        elif expected[0] == observed[0]:
            return 2  # correct
        else:
            return 3  # error
    elif is_missing and has_ref:
        return 0  # missing value
    return -1  # is_imputed = false and is not missing (non imputed value)

def compare_marker_between_two_lines(expected_genotypes, observed_genotypes, stats):
    samples_num = len(stats)
    assert len(expected_genotypes) == samples_num
    assert len(observed_genotypes) == samples_num
    for i in range(samples_num):
        r = compare_snp_marker_in_sample(expected_genotypes[i], observed_genotypes[i])
        if r != -1:
            stats[i][r] += 1
        if int(expected_genotypes[i][0]) < 3:
            stats[i][4] += 1

def update_lines(ref_line, test_line, expected_file, observed_file, stats):
    if not test_line:
        handle_missing_markers(ref_line, stats)
        ref_line = expected_file.readline().rstrip()
    elif not ref_line:
        handle_new_markers(test_line, stats)
        test_line = observed_file.readline().rstrip()
    else:
        assert ref_line and test_line
        ref_parts = ref_line.split('\t')
        test_parts = test_line.split('\t')
        if ref_parts[0] < test_parts[0]:
            handle_missing_markers(ref_line, stats)
            ref_line = expected_file.readline().rstrip()
        elif ref_parts[0] > test_parts[0]:
            handle_new_markers(test_line, stats)
            test_line = observed_file.readline().rstrip()
        else:
            assert ref_parts[0] == test_parts[0]
            compare_marker_between_two_lines(ref_parts[1:], test_parts[1:], stats)
            ref_line = expected_file.readline().rstrip()
            test_line = observed_file.readline().rstrip()
    return [ref_line, test_line]

def compare_files(expected, observed):
    expected_file = open(expected, "r")
    observed_file = open(observed, "r")
    ref_line = expected_file.readline().rstrip()
    test_line = observed_file.readline().rstrip()
    ref_parts = ref_line.split('\t')
    test_parts = test_line.split('\t')
    assert len(ref_parts) == len(test_parts)
    samples_num = len(ref_parts) - 1
    stats = [[0 for x in range(5)] for y in range(samples_num)]
    while ref_line or test_line:
        [ref_line, test_line] = update_lines(ref_line, test_line, expected_file, observed_file, stats)
    observed_file.close()
    expected_file.close()
    imputed_count = [(s[1]+s[2]+s[3]) for s in stats]
    missing_percentage = [(100.0*s[0])/s[4] for s in stats]
    error_percentage = [(100.0 * s[3]) / (s[1]+s[2]+s[3]) for s in stats]
    imputed_count_mean = np.mean(imputed_count)
    imputed_count_std = np.std(imputed_count)
    missing_percentage_mean = np.mean(missing_percentage)
    missing_percentage_std = np.std(missing_percentage)
    error_percentage_mean = np.mean(error_percentage)
    error_percentage_std = np.std(error_percentage)
    return [imputed_count_mean, imputed_count_std, missing_percentage_mean, missing_percentage_std, error_percentage_mean, error_percentage_std]


expected=sys.argv[1]
observed=sys.argv[2]
[imputed_count_mean, imputed_count_std, missing_percentage_mean, missing_percentage_std, error_percentage_mean,
 error_percentage_std] = compare_files(expected, observed)
str2 = "{0:.2f}".format(imputed_count_mean)
str3 = "{0:.2f}".format(imputed_count_std)
str4 = "{0:.2f}".format(missing_percentage_mean)
str5 = "{0:.2f}".format(missing_percentage_std)
str6 = "{0:.2f}".format(error_percentage_mean)
str7 = "{0:.2f}".format(error_percentage_std)
print "{}\t{}\t{}\t{}\t{}\t{}".format(str2, str3, str4, str5, str6, str7)

