#!/usr/bin/python
import sys,os
sys.path.append('/prodslow/testing/ariel/genotyper/code')
import GenotyperUtils

GenotyperUtils.snpToGenomagicInput(sys.argv[1],sys.argv[2])
