#!/usr/bin/python
import sys,os
sys.path.append('/prodslow/testing/ariel/genotyper/code')
import GenotyperUtils

GenotyperUtils.mergeSimilarityVCF(sys.argv[1],sys.argv[2])
#'simulations_of_1200_F1_1200.vcf','simulations_of_1200_F1_1200_short.vcf')

