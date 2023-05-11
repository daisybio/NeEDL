#!/bin/bash

ALL="../res/MACOED/"

python aggregate_results.py --csv-directories ${ALL} --save-as ../res/MACOED/ALL.csv --recursive
python aggregate_results.py --csv-directories ${ALL}ME/ --save-as ../res/MACOED/ME.csv
python aggregate_results.py --csv-directories ${ALL}NME/ --save-as ../res/MACOED/NME.csv