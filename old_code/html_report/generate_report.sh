#!/bin/bash
my_dir=$(dirname "$0")
my_data_version=$1
my_notebook=$my_dir/data_version_report.ipynb
export DATA_VERSION=$my_data_version
jupyter nbconvert $my_notebook --TemplateExporter.exclude_input=True --ExecutePreprocessor.enabled=True --ExecutePreprocessor.allow_errors=True --ExecutePreprocessor.timeout=180 --output-dir=./ --to html



