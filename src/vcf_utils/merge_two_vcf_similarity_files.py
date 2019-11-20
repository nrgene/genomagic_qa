import sys
import vcf_pair_iterator

def my_simple_func(line1, line2):
    if line1[:27]==line2[:27]:
        print(line1[:27])

progeny_vcf_file = sys.argv[1]
parents_vcf_file = sys.argv[2]
first_file = open(progeny_vcf_file, 'r')
second_file = open(parents_vcf_file, 'r')
line1 = vcf_pair_iterator.get_first_non_header_line(first_file)
line2 = vcf_pair_iterator.get_first_non_header_line(second_file)
vcf_pair_iterator.iterate_two_vcf(first_file, second_file, line1, line2, my_simple_func)
first_file.close()
second_file.close()