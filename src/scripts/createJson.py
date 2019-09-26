import sys


def create_sample_string(sample_name, path):
    return "\"{}\":[\"{}/{}_R1.fastq\",\"{}/{}_R2.fastq\"]".format(sample_name, path, sample_name, path, sample_name)


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


def create_population_of_parents(parent1, parent2, progeny_name, progeny_count, path, prefix_list):
    parentalLines = "\"parentalLines\":[\"{}\",\"{}\"]".format(parent1, parent2)
    sampleToFiles = "\"sampleToFiles\":{}\n{}\n{}".format("{", create_all_progeny(progeny_name, progeny_count, path, prefix_list),"}")
    return "{}{},{}{}".format("{", parentalLines, sampleToFiles, "}")


def create_all_populations(path, prefix_list):
    population_string1 = create_population_of_parents("b73v4__ver100", "lh82__ver110", "b73_lh82", 16, path, prefix_list)
    population_string2 = create_population_of_parents("b73v4__ver100", "ph207__ver100", "b73_ph207", 16, path,
                                                      prefix_list)
    population_string3 = create_population_of_parents("lh82__ver110", "ph207__ver100", "lh82_ph207", 20, path,
                                                      prefix_list)
    population_string4 = create_population_of_parents("b73v4__ver100", "phg39__ver110", "b73_phg39", 40, path,
                                                      prefix_list)
    return "\"populations\":[\n{},\n{},\n{},\n{}]".format(population_string1, population_string2, population_string3, population_string4)


def print_full_json(path, prefix_list, source, target, pedigree):
    sourceDataVersion = "\"sourceDataVersion\":\"{}\"".format(source)
    targetDataVersion = "\"targetDataVersion\":\"{}\"".format(target)
    pedigreeFileId = "\"pedigreeFileId\":\"{}\"".format(pedigree)
    populations_string = create_all_populations(path, prefix_list)
    return "{}\n{},\n{},\n{},\n{}\n{}".format("{",sourceDataVersion,targetDataVersion, pedigreeFileId,populations_string, "}")


prefix_list = []
for i in range(20):
    prefix_list.append("c{}".format(i+1))
path = sys.argv[4]
source_data_version = sys.argv[1]
target_data_version = sys.argv[2]
pedigree_file = sys.argv[3]
json_string = print_full_json(path, prefix_list, source_data_version, target_data_version, pedigree_file)
print json_string