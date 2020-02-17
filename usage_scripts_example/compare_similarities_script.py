import api_server.requests as ar
import analysis.similarities_comparison as sc

api_server = 'api-dev.nrgene.local:8080'
dv1='maize_benchmark_test_fix_mkrs_919_03'
dv2 = 'maize_benchmark_only_arg_wgs'
dv3 = 'dm_gm_public_maize_232'
samples_wgs=['b73v4__ver100', 'cml247__ver100', 'ep1_v2__ver100', 'w22__ver100', 'ki3__ver110', 'f7_v2__ver100', 'mo17__ver100']
samples_snp = ['b73', 'cml247', 'ep1', 'w22', 'ki3', 'f7', 'mo17']
n = len(samples_wgs)
my_dict = {}
for i in range(n):
    my_dict[samples_snp[i]] = samples_wgs[i]
hap_sim1 = ar.get_raw_similarities_between_multiple_sampels(api_server, dv1, samples_wgs, "1,2")
hap_sim2 = ar.get_raw_similarities_between_multiple_sampels(api_server, dv2, samples_wgs, "1,2")
snp_sim = ar.get_raw_similarities_between_multiple_sampels(api_server, dv3, samples_snp, "1,2").replace(my_dict)
[fn1, fp1, intersect1] = sc.compute_similarity_match_for_multiple_sampels_in_df(hap_sim1, snp_sim, samples_wgs)
[fn2, fp2, intersect2] = sc.compute_similarity_match_for_multiple_sampels_in_df(hap_sim2, snp_sim, samples_wgs)
print("{} {} {}".format(fn1/1000000, fp1/1000000, intersect1/1000000))
print("{} {} {}".format(fn2/1000000, fp2/1000000, intersect2/1000000))