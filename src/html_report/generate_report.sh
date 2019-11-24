#!/bin/bash
my_dir=$(dirname "$0")
jupyter nbconvert $my_dir/test.ipynb --output-dir=.


