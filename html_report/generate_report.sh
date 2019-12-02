#!/bin/bash
my_dir=$(dirname "$0")
my_data_version=$1
my_rs_host=$2
my_notebook=$my_dir/data_version_report.ipynb
export DATA_VERSION=$my_data_version
export REDSHIFT_HOST=$my_rs_host
jupyter nbconvert $my_notebook --TemplateExporter.exclude_input=True --ExecutePreprocessor.enabled=True --ExecutePreprocessor.allow_errors=True --ExecutePreprocessor.timeout=180 --output-dir=./ --to html



