import urllib3
import pandas as pd
import similarities_comparison

def api_call_request(url):
    http = urllib3.PoolManager()
    r = http.request('GET', url)
    return str(r.data).replace("\\t", "\t").replace("\\n", "\n")


def snp_marker_request(api_server, data_version, samples, locations):
    samples_query = samples.replace(",", "+OR+")
    url = "http://{}:8080/genomagic-api/v1/sync-job/SNP_MARKERS_COLOR_VCF.vcf?dataVersion={}&samplesQuery={}" \
          "&locations={}&outputSamples={}".format(api_server, data_version, samples_query, locations, samples)
    my_data = api_call_request(url)
    print(my_data)


def genetic_distance_request(api_server, data_version, samples, locations):
    url = "http://{}:8080/genomagic-api/v1/sync-job/HAPLOTYPES_SIMILARITY_GENETIC_DISTANCE.tsv?dataVersion={}&samples={}" \
          "&locations={}".format(api_server, data_version, samples, locations)
    return api_call_request(url)


def get_genetic_distance_values_from_api(api_server, data_version, samples, locations):
    data = genetic_distance_request(api_server, data_version, samples, locations)
    lines = data[2:-2].split("\n")
    sample1 = [x.split("\t")[0] for x in lines]
    sample2 = [x.split("\t")[1] for x in lines]
    val = [float(x.split("\t")[2]) for x in lines]
    df = pd.DataFrame(sample1, columns=['sample 1'])
    df['sample 2'] = sample2
    df['genetic distance'] = val
    return df


def get_raw_similarities_between_two_sampels(api_server, data_version, sample1, sample2, locations):
    samples = "{},{}".format(sample1, sample2)
    df = get_raw_similarities_between_multiple_sampels(api_server, data_version, samples, locations)
    return df.loc[:,['chr', 'start', 'end']]


def get_raw_similarities_between_multiple_sampels(api_server, data_version, samples, locations):
    url = "http://{}/genomagic-api/v1/sync-job/HAPLOTYPES_SIMILARITY_RAW.tsv?dataVersion={}&samples={}" \
          "&locations={}".format(api_server, data_version, ",".join(samples), locations)
    print(url)
    my_data = api_call_request(url)
    lines = my_data[2:-2].split("\n")
    data = [x.split('\t') for x in lines]
    df = pd.DataFrame(data)
    df.columns = ["s1", "s2", "chr", "start", "end", "score", "ibd"]
    ind = (df["ibd"] == "true") & (df["start"] != df["end"])
    hap_sim = df.loc[ind,['s1','s2','chr', 'start', 'end']]
    hap_sim['chr'] = hap_sim.chr.astype(int)
    hap_sim['start'] = hap_sim.start.astype(int)
    hap_sim['end'] = hap_sim.end.astype(int)
    return hap_sim


def get_chromosome_bounds(api_server, data_version, sample):
    url = 'http://{}/genomagic-api/v1/chromosomes-bounds?dataVersion={}&sampleId={}&responseType=tsv'.format(api_server, data_version, sample)
    my_data = api_call_request(url)
    lines = my_data[2:-2].split("\n")
    data = [x.split('\t') for x in lines]
    df = pd.DataFrame(data)
    df.columns = ["chr", "start", "end"]
    return df.astype(int)


def parse_locations(locations_string, chromosome_bounds):
    locations = locations_string.split(',')
    int_locations = []
    for location in locations:
        loc_sep = location.split(':')
        my_chr = int(loc_sep[0])
        if len(loc_sep) ==1:
            my_chr_bounds = chromosome_bounds.loc[chromosome_bounds['chr'] == my_chr]
            assert my_chr_bounds.shape == (1, 3)
            range_start = int(my_chr_bounds.start)
            range_end = int(my_chr_bounds.end)
        else:
            assert len(loc_sep) == 2
            my_range = loc_sep[1].split('-')
            range_start = int(my_range[0])
            range_end = int(my_range[1])
        int_locations.append([my_chr, range_start, range_end])
    return int_locations


def compute_genetic_distance_between_two_samples(api_server, data_version, samples, locations_string, pivot_sample):
    chromosome_bounds = []
    all_locations_have_start_and_end = locations_string.count("-") == locations_string.count(",") + 1
    if not all_locations_have_start_and_end:
        chromosome_bounds = get_chromosome_bounds(api_server, data_version, pivot_sample)

    loc = parse_locations(locations_string, chromosome_bounds)
    hap_sim = get_raw_similarities_between_multiple_sampels(api_server, data_version, samples,
                                                                             locations_string)

    all_merged_clipped_list = []
    total_len = 0
    similarity_len = 0

    for i in range(len(loc)):
        total_len += loc[i][2] - loc[i][1]
        ind = (hap_sim.chr == loc[i][0]) & (hap_sim.start <= loc[i][2]) & (hap_sim.end >= loc[i][1])
        temp_df = hap_sim.loc[ind, ['start', 'end']]
        merged_list = similarities_comparison.merge_list_elements(temp_df.values)
        merged_clipped_list = [[loc[i][0], max(x[0], loc[i][1]), min(x[1], loc[i][2])] for x in merged_list]
        all_merged_clipped_list = all_merged_clipped_list + merged_clipped_list

    for similarity in all_merged_clipped_list:
        similarity_len += similarity[2] - similarity[1]

    return round((1.0 * similarity_len) / total_len, 2)









