import urllib3
import pandas as pd


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
    return df.loc[ind,['s1','s2','chr', 'start', 'end']]

    # df = df.transpose()
    #df.columns = ["s1", "s2", "chr", "start", "end", "score", "ibd"]
    #ind = (df["ibd"] == "true") & (df["start"] != df["end"])
    #df = df.loc[ind, ['chr', 'start', 'end']].astype(int).drop_duplicates()
    #return df




