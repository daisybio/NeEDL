#!/bin/bash

ROOT="../res/EpiGEN"
DICHO="dichotomous"
CATEG="categorical"
QUANT="quantitative"
SIZE2="2_disease_snps"
SIZE3="3_disease_snps"
SIZE4="4_disease_snps"
DOM="joint-dominant"
REC="joint-recessive"
MULT="multiplicative"
EXP="exponential"

#for TYPE in ${DICHO} ${QUANT}
#do
#	python aggregate_results.py --csv-directories ${ROOT}/${TYPE}/ --save-as ${ROOT}/${TYPE}.csv --recursive
#	for SIZE in ${SIZE2} ${SIZE3} ${SIZE4}
#	do
#		python aggregate_results.py --csv-directories ${ROOT}/${TYPE}/${SIZE}/ --save-as ${ROOT}/${TYPE}_${SIZE}.csv --recursive
#	done
#	for MODEL in ${DOM} ${REC} ${MULT} ${EXP}
#	do
#		DIRS="${ROOT}/${TYPE}/${SIZE2}/${MODEL}/ ${ROOT}/${TYPE}/${SIZE3}/${MODEL}/ ${ROOT}/${TYPE}/${SIZE4}/${MODEL}/"
#		python aggregate_results.py --csv-directories ${DIRS} --save-as ${ROOT}/${TYPE}_${MODEL}.csv
#	done
#done

for TYPE in ${CATEG}
do
	for MODEL in "model_1" "model_2"
	do
		DIRS="${ROOT}/${TYPE}/${SIZE2}/${MODEL}/ ${ROOT}/${TYPE}/${SIZE3}/${MODEL}/ ${ROOT}/${TYPE}/${SIZE4}/${MODEL}/"
		python aggregate_results.py --csv-directories ${DIRS} --save-as ${ROOT}/${TYPE}_${MODEL}.csv
	done
done