def create_sample_string(sample_name, file1, file2):
    return "\"{}\":[\"{}\",\"{}\"]".format(sample_name, file1, file2)


def create_multiple_samples_string(arr):
    my_str = ''
    for x in arr:
        assert len(x) == 3
        sample_string = create_sample_string(x[0], x[1], x[2])
        my_str = '{},\n{}'.format(my_str, sample_string)
    return my_str[2:]


def create_sample_to_files_string(progeny_files):
    multiple_samples_string = create_multiple_samples_string(progeny_files)
    return "\"sampleToFiles\":{}\n{}\n{}".format("{", multiple_samples_string, "}")


def create_population_of_parents(parent1, parent2, progeny_files):
    parentalLines = "\"parentalLines\":[\"{}\",\"{}\"]".format(parent1, parent2)
    sampleToFiles = create_sample_to_files_string(progeny_files)
    return "{}{},{}{}".format("{", parentalLines, sampleToFiles, "}")


def create_populations_string_without_pedigree_file(parents_list, progeny_files_list):
    n = len(parents_list)
    assert len(progeny_files_list) == n
    all_parents_list = []
    for i in range(n):
        parents_pair = parents_list[i]
        assert len(parents_pair) == 2
        parents_str = create_population_of_parents(parents_pair[0], parents_pair[1], progeny_files_list[i])
        all_parents_list.append(parents_str)
    all_parents_string = ",\n".join(all_parents_list)
    return "\"populations\":[\n{}]".format(all_parents_string)


def create_population_string(pedigreeFileId, parents_list, progeny_files_list):
    populations_string_without_pedigree_file = create_populations_string_without_pedigree_file(parents_list, progeny_files_list)
    return "\"pedigreeFileId\":\"{}\",\n{}".format(pedigreeFileId, populations_string_without_pedigree_file)


def get_full_json_str(prev_dv, new_dv, pop_str):
    dv_str = "\"sourceDataVersion\":\"{}\",\n\"sourceDataVersion\":\"{}\",".format(prev_dv, new_dv)
    return "{}\n{}\n{}\n{}".format("{", dv_str, pop_str, "}")


def create_parental_json(prev_dv, new_dv, pedigreeFileId, parents_list, progeny_files_list):
    pop_str = create_population_string(pedigreeFileId, parents_list, progeny_files_list)
    return get_full_json_str(prev_dv, new_dv, pop_str)


def create_diversity_json(prev_dv, new_dv, progeny_files):
    pop_str = create_sample_to_files_string(progeny_files)
    return get_full_json_str(prev_dv, new_dv, pop_str)


def create_all_prefix_for_sample(sample_name, prefix_list, path):
    arr = []
    for prefix in prefix_list:
        new_sample_name = "{}_{}".format(prefix, sample_name)
        arr.append(create_sample_string(new_sample_name, path))
    return ",\n".join(arr)


def create_all_progeny(sample_name, number_of_progeny, path, prefix_list):
    arr = []
    for i in range(number_of_progeny):
        arr.append(create_all_prefix_for_sample("{}_{}".format(sample_name, i+1), prefix_list, path))
    return ",\n".join(arr)





def create_all_populations(path, prefix_list):
    population_string1 = create_population_of_parents("b73v4__ver100", "lh82__ver110", "b73_lh82", 16, path, prefix_list)
    population_string2 = create_population_of_parents("b73v4__ver100", "ph207__ver100", "b73_ph207", 16, path,
                                                      prefix_list)
    population_string3 = create_population_of_parents("lh82__ver110", "ph207__ver100", "lh82_ph207", 20, path,
                                                      prefix_list)
    population_string4 = create_population_of_parents("b73v4__ver100", "phg39__ver110", "b73_phg39", 40, path,
                                                      prefix_list)
    #return "\"populations\":[\n{},\n{},\n{},\n{}]".format(population_string1, population_string2, population_string3, population_string4)
    return "\"populations\":[\n{},\n{},\n{}]".format(population_string1, population_string3, population_string4) 

def print_full_json(path, prefix_list, source, target, pedigree):
    sourceDataVersion = "\"sourceDataVersion\":\"{}\"".format(source)
    targetDataVersion = "\"targetDataVersion\":\"{}\"".format(target)
    pedigreeFileId = "\"pedigreeFileId\":\"{}\"".format(pedigree)
    populations_string = create_all_populations(path, prefix_list)
    return "{}\n{},\n{},\n{},\n{}\n{}".format("{",sourceDataVersion,targetDataVersion, pedigreeFileId,populations_string, "}")

def print_full_json_diversity_analysis(path, prefix_list, source, target):
    sourceDataVersion = "\"sourceDataVersion\":\"{}\"".format(source)
    targetDataVersion = "\"targetDataVersion\":\"{}\"".format(target)
    population1=create_all_progeny("b73_lh82".upper(), 16, path, prefix_list)
    population2=create_all_progeny("b73_ph207".upper(), 16, path, prefix_list)
    population3=create_all_progeny("lh82_ph207".upper(), 20, path, prefix_list)
    population4=create_all_progeny("b73_phg39".upper(), 40, path, prefix_list)
    all_files_string="{},\n{},\n{}".format(population1,population3,population4)    
    sampleToFiles = "\"sampleToFiles\":{}\n{}\n{}".format("{", all_files_string,"}")
    #populations_string = create_all_populations(path, prefix_list)
    return "{}\n{},\n{},\n{}\n{}".format("{",sourceDataVersion,targetDataVersion,sampleToFiles,"}")


#s = create_string_for_single_sample("arr", "s3://file1.fastq", "s3://file1.fastq")
arr1 = []
arr1.append(["arr", "s3://file1.fastq", "s3://file2.fastq"])
arr1.append(["das", "s3://file3.fastq", "s3://file4.fastq"])
arr1.append(["ghh", "s3://file5.fastq", "s3://file6.fastq"])


arr2 = []
arr2.append(["asc", "s3://file1.fastq", "s3://file2.fastq"])
arr2.append(["A#$", "s3://file3.fastq", "s3://file4.fastq"])
arr2.append(["YYY", "s3://file5.fastq", "s3://file6.fastq"])


s3 = create_diversity_json("public_fake_v0_13", "public_fake_v0_14", arr1)



#s3 = create_population_of_parents("aviad", "amitay", arr1)
print(s3)
#print(s1)
#print(s2)

#s=[]
#print(s)
#import pandas as pd
#df = pd.read_csv('gbs_files_with_pedigree', sep='\t')
#df['p1'] = df[['parent1','parent2']].min(axis=1)
#df['p2'] = df[['parent1','parent2']].max(axis=1)

#pd.read_csv(fpath, sep='\t')


#DataFrame.from_csv('c:/~/trainSetRel3.txt', sep='\t')
#prefix_list = []
#for i in range(20):
#    prefix_list.append("c{}".format(i+1))
#path = sys.argv[4]
#source_data_version = sys.argv[1]
#target_data_version = sys.argv[2]
#pedigree_file = sys.argv[3]
#json_string = print_full_json(path, prefix_list, source_data_version, target_data_version, pedigree_file)
#json_string = print_full_json_diversity_analysis(path, prefix_list, source_data_version, target_data_version)
#print(json_string)
#print create_all_progeny("b73_lh82", 16, path, prefix_list)

