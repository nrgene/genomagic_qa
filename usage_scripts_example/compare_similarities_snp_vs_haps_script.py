import sys
import pandas as pd
genomagic_path = '/home/ariel/clients/genomagic_qa/'
sys.path.append(genomagic_path)
import redshift.redshift_queries as rq
import api_server.requests as api
import analysis.markers_comparison as mc
import analysis.similarities_comparison as sc

# hard coded parameters
api_server='api-dev.nrgene.local'
locations="1,2,3,4,5"
max_hap_num = 1000000
data_version = 'public_soy_v1_14'
host='rndlab-genomagic-redshift.cl6ox83ermwm.us-east-1.redshift.amazonaws.com'
samples = 'lee__ver100,pi483463__ver100,williams82__ver100,fc33243_anderson__ver100,pi180501_strain_no__18__ver100,pi240664_bilomi_no__3__ver100,pi266806c_no_4__ver100,pi322692_max_c_p1159a8__ver100,pi360957_karafuto_no__1__ver100,pi361087_medias_23__ver100,pi361093_novosadska_br__1__ver100,pi437265d_dobruzanca_d__ver100,pi438335_sao_196_c__ver100,pi468908__ver100,pi497964a_i_c__9461__ver100,pi507180_rikuu_21__ver100,pi518668_tn_4_86__ver100,pi548193_t201__ver100,pi548360_korean__ver100,pi548364_macoupin__ver100,pi548447_cherokee__ver100,pi548452_dixie__ver100,pi548490_tanner__ver100,pi548520_preston__ver100,pi548561_hodgson__ver100,pi548619_sparks__ver100,pi549041a_zyd_2709__ver100,pi559932_manokin__ver100,pi567426_bai_huang_dou__ver100,pi567558_liu_shi_ri_jin_huang_da_dou__ver100,pi567604a_xin_huang_dou__ver100,pi567788_bienville__ver100,pi578309_i_64__ver100,pi592523_glacier__ver100'
similarity_threshold = 0.85
max_major_allele_freq = 0.8
min_samples_presence = 2
min_p = 1
win_len = 20
trim_len = 0


"""
In this script we are extracting snp data from wgs samples, and then generate similarities based on these snps
These similarities are compared to haplotype based similarities from redshift
The assumption is that this is a way to optimize both hap similarity parameters and markers similarity parameteres
"""


# here we get a pair of samples , and return the similarities comparison between snp similarities and hap similarities
def compare_similarities_for_samples_pair(sample1, sample2, hap_similarities, curr_snp_markers, informative, min_p, win_len, trim_len):
    ind = (hap_similarities['s1'] == sample1) & (hap_similarities['s2'] == sample2)
    hap_sim_list = sc.merge_list_elements(list(hap_similarities.loc[ind, ['start', 'end']].values))
    df_test = pd.DataFrame({'chr':curr_snp_markers['chr'], 'pos': curr_snp_markers['pos'],
                            'informative': informative, 'sample1': curr_snp_markers[sample1],
                            'sample2': curr_snp_markers[sample2]})
    [snp_sim_list, match_pr, mismatch_pr] = mc.compute_similarities_between_samples(df_test, min_p, win_len, trim_len)
    return sc.compare_lists_inrtersect(hap_sim_list, snp_sim_list)


# here for a given chromosome, we compare the snp vs hap similarities over all sample pairs
def compare_all_similarities_by_chr(samples_list, chromosome, host, similarity_table, similarity_threshold, snp_markers, max_major_allele_freq, min_samples_presence, min_p, win_len, trim_len):
    query = 'select sample1, sample2, chromosome_id, start_position, end_position from {} where chromosome_id={} and similarity_score>={} order by start_position;'.format(similarity_table, chromosome, similarity_threshold)
    hap_similarities = rq.get_sql_query_as_data_frame(host, query, ['s1','s2','chr','start','end'])
    samples_count = len(samples_list)
    curr_snp_markers = snp_markers[snp_markers['chr']==chromosome].reset_index(drop=True)
    informative = mc.compute_informative_markers(snp_markers, max_major_allele_freq, min_samples_presence)
    hap_sim_total = 0
    snp_sim_total = 0
    shared_sim_total = 0
    for i in range(samples_count):
        for j in range(i + 1, samples_count):
            sample1 = samples_list[i]
            sample2 = samples_list[j]
            [only_hap_sim, only_snp_sim, hap_snp_shared_sim]= compare_similarities_for_samples_pair(sample1, sample2, hap_similarities, curr_snp_markers, informative, min_p, win_len, trim_len)
            hap_sim_total += only_hap_sim
            snp_sim_total += only_snp_sim
            shared_sim_total += hap_snp_shared_sim
    return [hap_sim_total, snp_sim_total, shared_sim_total]


# here for a given data version, we compare the similarities between all samples in all chromosomes
def compare_sim_all_chr(samples_list, chr_string, host, data_version, similarity_threshold, snp_markers,
                        max_major_allele_freq, min_samples_presence, min_p, win_len, trim_len):
    similarity_table = rq.get_table_name(host, data_version, 'HAPLOTYPES_SIMILARITY')
    only_hap = 0
    only_snp = 0
    both_methods = 0
    chr_list = chr_string.split(',')
    for c in chr_list:
        chromosome = int(c)
        [hap_sim_total, snp_sim_total, shared_sim_total] = compare_all_similarities_by_chr(samples_list, chromosome,
                                                                                           host, similarity_table,
                                                                                           similarity_threshold,
                                                                                           snp_markers,
                                                                                           max_major_allele_freq,
                                                                                           min_samples_presence, min_p,
                                                                                           win_len, trim_len)
        only_hap += hap_sim_total
        only_snp += snp_sim_total
        both_methods += shared_sim_total
    return [only_hap, only_snp, both_methods]


snp_markers = api.get_snp_markers_as_dataframe(api_server, data_version, samples, locations, max_hap_num)
samples_list = samples.split(',')

data_version = 'public_soy_v1_14'
[only_hap, only_snp, both] = compare_sim_all_chr(samples_list, locations, host, data_version, similarity_threshold, snp_markers,
                        max_major_allele_freq, min_samples_presence, min_p, win_len, trim_len)
total_len = only_hap + only_snp + both
print('maxHaplotypeFrequency =  1.0')
print('dv = {} , jacard = {}, only hap = {} only snp = {}'.format(data_version, both/total_len, only_hap/total_len, only_snp/total_len))

data_version = 'public_soy_test_remove_freq'
[only_hap, only_snp, both] = compare_sim_all_chr(samples_list, locations, host, data_version, similarity_threshold, snp_markers,
                        max_major_allele_freq, min_samples_presence, min_p, win_len, trim_len)
total_len = only_hap + only_snp + both
print('maxHaplotypeFrequency =  0.9')
print('dv = {} , jacard = {}, only hap = {} only snp = {}'.format(data_version, both/total_len, only_hap/total_len, only_snp/total_len))

data_version = 'public_soy_test_remove_freq_3'
[only_hap, only_snp, both] = compare_sim_all_chr(samples_list, locations, host, data_version, similarity_threshold, snp_markers,
                        max_major_allele_freq, min_samples_presence, min_p, win_len, trim_len)
total_len = only_hap + only_snp + both
print('maxHaplotypeFrequency =  0.8')
print('dv = {} , jacard = {}, only hap = {} only snp = {}'.format(data_version, both/total_len, only_hap/total_len, only_snp/total_len))

data_version = 'public_soy_test_remove_freq_2'
[only_hap, only_snp, both] = compare_sim_all_chr(samples_list, locations, host, data_version, similarity_threshold, snp_markers,
                        max_major_allele_freq, min_samples_presence, min_p, win_len, trim_len)
total_len = only_hap + only_snp + both
print('maxHaplotypeFrequency =  0.6')
print('dv = {} , jacard = {}, only hap = {} only snp = {}'.format(data_version, both/total_len, only_hap/total_len, only_snp/total_len))




