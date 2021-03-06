import urllib3
import pandas as pd
import sys
genomagic_qa_repo_path = '/'.join(sys.argv[0].split('/')[:-2])
sys.path.append(genomagic_qa_repo_path)
import vcf_utils.vcf_reader as vr


def register_data_version_command(host, dv_name):
    url = 'http://{}/genomagic-api/v1/data-versions/{}'.format(host, dv_name)
    return 'curl -X PUT \'{}\''.format(url)


def add_samples_group_command(samples_list, group_name, host, dv_name):
    samples_string = ",".join(["\"{}\"".format(s) for s in samples_list])
    sample_string_json = "{}\"sampleIds\" : [{}]{}".format("{", samples_string, "}")
    url = "http://{}/genomagic-api/v1/sample-groups/{}?dataVersion={}".format(host, group_name, dv_name)
    my_command = " curl -X PUT \'{}\' - H \'Content-Type: application/json\' -d \'{}\'".format(url, sample_string_json)
    return my_command


def api_call_request(url):
    http = urllib3.PoolManager()
    r = http.request('GET', url)
    if r.status == 202:
        print("request didn't finish yet")
        assert r.status == 200, r.data
    elif r.status != 200:
        assert r.status == 200, r.data
    return str(r.data).replace("\\t", "\t").replace("\\n", "\n")


def snp_marker_request(api_server, data_version, samples, locations):
    samples_query = samples.replace(",", "+OR+")
    url = "http://{}:8080/genomagic-api/v1/sync-job/SNP_MARKERS_COLOR_VCF.vcf?dataVersion={}&samplesQuery={}" \
          "&locations={}&outputSamples={}".format(api_server, data_version, samples_query, locations, samples)
    my_data = api_call_request(url)
    lines = my_data[2:-2].split("\n")
    return lines


def haplotype_marker_request(api_server, data_version, samples, locations):
    samples_query = samples.replace(",", "+OR+")
    url = "http://{}:8080/genomagic-api/v1/sync-job/HAPLOTYPE_MARKERS_COLOR_VCF.vcf?dataVersion={}&samplesQuery={}" \
          "&locations={}&outputSamples={}".format(api_server, data_version, samples_query, locations, samples)
    print(url)
    my_data = api_call_request(url)
    lines = my_data[2:-2].split("\n")
    return lines


def genetic_distance_request(api_server, data_version, samples, locations):
    url = "http://{}:8080/genomagic-api/v1/sync-job/HAPLOTYPES_SIMILARITY_GENETIC_DISTANCE.tsv?dataVersion={}&samples={}" \
          "&locations={}".format(api_server, data_version, samples, locations)
    return api_call_request(url)


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


def get_snp_markers_as_dataframe(api_server, data_version, samples, locations, max_haps_num):
    vcf_lines_list = snp_marker_request(api_server, data_version, samples, locations)
    return vr.read_snp_vcf_from_list(vcf_lines_list, max_haps_num)


def get_haplotype_markers_as_dataframe(api_server, data_version, samples, locations, max_haps_num):
    vcf_lines_list = haplotype_marker_request(api_server, data_version, samples, locations)
    return vr.read_hap_vcf_from_list(vcf_lines_list, max_haps_num)



