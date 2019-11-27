from urllib2 import Request, urlopen
from operator import itemgetter
import time


def test_imputation_hap(comparison_list, my_api_server, my_data_version, locations_string,
                        false_positive_threshold, false_negative_threshold):
    imputed_samples = [x[0] for x in comparison_list]
    original_samples = [x[1] for x in comparison_list]
    markers_matching_stats = compare_hap_markers_imputation(my_api_server, my_data_version,
                                                                             locations_string, original_samples,
                                                                             imputed_samples)
    #print markers_matching_stats
    assert len(original_samples) == len(imputed_samples)
    no_false_positive = test_value(markers_matching_stats, 1, false_positive_threshold, original_samples,
                                   imputed_samples, "false_positive")
    no_false_negative = test_value(markers_matching_stats, 2, false_negative_threshold, original_samples,
                                   imputed_samples, "false_negative")
    if no_false_positive and no_false_negative:
        print 'TESTS PASSED'


def test_imputation_snp(comparison_list, my_api_server, my_data_version, locations_string, mismatch_percent_threshold,
                        missing_percent_threshold, new_percent_threshold):
    imputed_samples = [x[0] for x in comparison_list]
    original_samples = [x[1] for x in comparison_list]
    markers_matching_stats = compare_snp_markers_imputation(my_api_server, my_data_version, locations_string,
                                                            original_samples, imputed_samples)
    #print markers_matching_stats
    assert len(original_samples) == len(imputed_samples)
    no_mismatch = test_value(markers_matching_stats, 0, mismatch_percent_threshold, original_samples, imputed_samples, "mismatch")
    no_new = test_value(markers_matching_stats, 1, new_percent_threshold, original_samples, imputed_samples, "new")
    no_missing = test_value(markers_matching_stats, 2, missing_percent_threshold, original_samples, imputed_samples, "missing")
    if no_mismatch and no_new and no_missing:
        print 'TESTS PASSED'


def test_value(markers_matching_stats, index, threshold, original_samples, imputed_samples, test_name):
    samples_count = len(original_samples)
    test_success = True
    for i in range(samples_count):
        total_markers = sum(markers_matching_stats[i]) * 1.0
        observed_value = round(markers_matching_stats[i][index] / total_markers * 100, 2)
        if observed_value > threshold:
            test_success = False
            print 'FAILED: {} imputation vs {} source {} = {}%' \
                .format(imputed_samples[i], original_samples[i], test_name, observed_value)
    return test_success






def compare_hap_markers_imputation(my_api_server, my_data_version, locations_string, original_samples, imputed_samples):
    start_time = time.time()
    print "sending api requests - haplotypes"
    markers_file = marker_request(my_api_server, my_data_version, ','.join(original_samples), locations_string,
                                  "HAPLOTYPE")
    imputation_file = imputation_request(my_api_server, my_data_version, ','.join(imputed_samples), locations_string,
                                         'haplotype')
    end_time = time.time()
    print "api request finished successfully after {} seconds".format(round(end_time - start_time,1))
    markers_matching_stats = compare_markers_file(markers_file, imputation_file, 4, compare_hap_marker)
    return markers_matching_stats


def compare_snp_markers_imputation(my_api_server, my_data_version, locations_string, original_samples, imputed_samples):
    start_time = time.time()
    print "sending api request - snp"
    markers_file = marker_request(my_api_server, my_data_version, ','.join(original_samples), locations_string, "SNP")
    imputation_file = imputation_request(my_api_server, my_data_version, ','.join(imputed_samples), locations_string,
                                         'snp')
    end_time = time.time()
    print "api request finished successfully after {} seconds".format(round(end_time - start_time,1))
    markers_matching_stats = compare_markers_file(markers_file, imputation_file, 2, compare_snp_marker)
    return markers_matching_stats


def compare_markers_file(file1, file2, id_index, compare_method):
    markers_list1 = read_markers_file(file1, id_index)
    markers_list2 = read_markers_file(file2, id_index)
    markers1_count = len(markers_list1)
    markers2_count = len(markers_list2)
    assert len(markers_list1[0]) == len(markers_list2[0])
    samples_num = len(markers_list1[0])-1
    # hap stats[0]: missing in both, stats[1] only in imp , stats[2] only in hap , stats[3] both
    # snp stats[0]: mismatch, stats[1] only in imp , stats[2] only in hap , stats[3] match
    stats = [[0 for x in range(4)] for y in range(samples_num)]
    marker1_index = 0
    marker2_index = 0
    while marker1_index < markers1_count or marker2_index < markers2_count:
        if marker1_index < markers1_count:
            marker1_id = markers_list1[marker1_index][0]
        if marker2_index < markers2_count:
            marker2_id = markers_list2[marker2_index][0]
        if marker1_id > marker2_id or marker1_index >= markers1_count:  # missing first marker
            for i in range(samples_num):
                stats[i][1] += 1
            marker2_index += 1
        elif marker2_id > marker1_id or marker2_index >= markers2_count:  # missing second marker
            for i in range(samples_num):
                stats[i][2] += 1
            marker1_index += 1
        else:  # same marker
            assert marker1_id == marker2_id
            assert marker1_index < markers1_count and marker2_index < markers2_count
            compare_method(markers_list1[marker1_index], markers_list2[marker2_index], stats)
            marker1_index += 1
            marker2_index += 1
    return stats


def read_markers_file(markers_file, id_column):
    line = markers_file.readline().rstrip()
    vcf_content = []
    while line[0] == "#":
        line = markers_file.readline().rstrip()
    while line:
        parts = line.split('\t')
        vcf_content.append(parts[id_column:id_column+1]+parts[9:])
        line = markers_file.readline().rstrip()
    return sorted(vcf_content, key=itemgetter(0))


def marker_request(api_server, data_version, samples, locations, job):
    samples_query = samples.replace(",", "+OR+")
    url = "http://{}:8080/genomagic-api/v1/sync-job/{}_MARKERS_COLOR_VCF.vcf?dataVersion={}&samplesQuery={}" \
          "&locations={}&outputSamples={}".format(api_server, job, data_version, samples_query, locations, samples)
    print url
    result = urlopen(Request(url))
    return result


def imputation_request(api_server, data_version, samples, locations, job):
    url = "http://{}:8080/genomagic-api/v1/sync-job/MARKER_IMPUTATION_COLOR_VCF.vcf?dataVersion={}&markerType={}" \
          "&samples={}&locations={}".format(api_server, data_version, job, samples, locations)
    print url
    result = urlopen(Request(url))
    return result


def compare_snp_marker(marker1, marker2, stats):
    assert len(marker1) == len(marker2), "{} samples in 1 {} samples in 2 - should be the same"\
        .format(len(marker1)-1, len(marker2)-1)
    assert marker1[0] == marker2[0]
    samples_num = len(marker1)-1
    for i in range(samples_num):
        val1 = int(marker1[i+1][0])
        val2 = int(marker2[i+1][0])
        if val2 == val1 <= 2:  # we have a match
            stats[i][3] += 1
        elif val2 <= 2 and val1 <= 2:  # we have a mismatch
            assert val2 != val1
            stats[i][0] += 1
        elif val1 > 2:  # missing first marker
            stats[i][1] += 1
        elif val2 > 2:  # missing second marker
            stats[i][2] += 1


def compare_hap_marker(marker1, marker2, stats):
    assert len(marker1) == len(marker2), "{} samples in 1 {} samples in 2 - should be the same"\
        .format(len(marker1)-1, len(marker2)-1)
    assert marker1[0] == marker2[0]
    samples_num = len(marker1)-1
    for i in range(samples_num):
        if marker1[i+1][0] == '1' and marker2[i+1][0] == '1':
            stats[i][0] += 1
        elif marker1[i+1][0] == '1' and marker2[i+1][0] == '2':
            stats[i][1] += 1
        elif marker1[i+1][0] == '2' and marker2[i+1][0] == '1':
            stats[i][2] += 1
        else:
            assert marker1[i+1][0] == '2' and marker2[i+1][0] == '2'
            stats[i][3] += 1



