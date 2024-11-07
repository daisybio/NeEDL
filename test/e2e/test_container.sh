#!/bin/bash

# This scripts tests the functionality of the docker container
# Please provide the docker image name as first parameter to this script.


# ---- setup

# abort on first error
set -e

# create temp output folder
mkdir -p ./test_out
out_dir=$(realpath ./test_out)

echo "Temp output directory: $out_dir"
echo "Working directory: $(pwd)"

if [ -z "$2" ]; then
    echo "Using native platform."
    py_platform_str=""
else
    echo "Using platform: $2"
    py_platform_str="--docker-platform $2"
fi

# ---- tests

# check if binaries exist and are executable
python ./run/NeEDL.py --container-image-name "$1" --docker-no-pulling $py_platform_str --help
python ./run/epiJSON.py --container-image-name "$1" --docker-no-pulling $py_platform_str --help
python ./run/calculate_scores.py --container-image-name "$1" --docker-no-pulling $py_platform_str --help
python ./run/convert_to_binary.py --container-image-name "$1" --docker-no-pulling $py_platform_str --help

# check if it also works if we explicitly specify docker as container platform
python ./run/NeEDL.py --docker --container-image-name "$1" --docker-no-pulling $py_platform_str --help
python ./run/epiJSON.py --docker --container-image-name "$1" --docker-no-pulling $py_platform_str --help
python ./run/calculate_scores.py --docker --container-image-name "$1" --docker-no-pulling $py_platform_str --help
python ./run/convert_to_binary.py --docker --container-image-name "$1" --docker-no-pulling $py_platform_str --help

# realtime_scores is a special case --> no python launcher script exists for it
# run containers as current user
user_id=$(id -u)
group_id=$(id -g)
docker run --user $user_id:$group_id "$1" /NeEDL/test/model/bin/realtime_scores --help



# select small dummy dataset
dummy_dataset=$(realpath ./data/e2e_tests/dummy_dataset)

# test NeEDL quick start example with dummy dataset
python ./run/NeEDL.py \
    --container-image-name "$1" \
    --docker-no-pulling \
    $py_platform_str \
    "--num-threads" "1" \
    "--output-directory" "$out_dir" \
    "--input-format" "JSON_EPIGEN" \
    "--input-path" "$dummy_dataset.json" \
    "--phenotype" "DICHOTOMOUS" \
    "--snp-annotate-dbSNP" \
    "--network-BIOGRID" \
    "--ms-seeding-routine" "RANDOM_CONNECTED" \
    "--ms-rc-start-seeds" "5" \
    "--ms-model" "PENETRANCE_NLL"


# check that eqtl_mapping works
python ./run/NeEDL.py \
    --container-image-name "$1" \
    --docker-no-pulling \
    $py_platform_str \
    "--num-threads" "1" \
    "--output-directory" "$out_dir" \
    "--input-format" "JSON_EPIGEN" \
    "--input-path" "$dummy_dataset.json" \
    "--phenotype" "DICHOTOMOUS" \
    "--snp-annotate-eQTL" \
    "--network-BIOGRID" \
    "--ms-seeding-routine" "RANDOM_CONNECTED" \
    "--ms-rc-start-seeds" "5" \
    "--ms-model" "PENETRANCE_NLL"


# test that custom annotation files and custom networks work
annotation_file=$(realpath ./data/e2e_tests/annotations.csv)
network_file=$(realpath ./data/e2e_tests/network.csv)

# currently this creates an empty network --> checks only if annotation and network file is processed correctly
# Maybe replace with something that actually creates a network from the custom data?
python ./run/NeEDL.py \
    --container-image-name "$1" \
    --docker-no-pulling \
    $py_platform_str \
    "--num-threads" "1" \
    "--output-directory" "$out_dir" \
    "--input-format" "JSON_EPIGEN" \
    "--input-path" "$dummy_dataset.json" \
    "--phenotype" "DICHOTOMOUS" \
    "--snp-annotate" "$annotation_file|yes|SNP|gene|\t|-1|-1" \
    "--network" "dummy_network|$network_file|yes|col1|col2|\t|-1|-1" \
    "--ms-seeding-routine" "RANDOM_CONNECTED" \
    "--ms-rc-start-seeds" "5" \
    "--ms-model" "PENETRANCE_NLL"
 

# epiJSON convert json --> all formats
epijson_test_dir_1=$out_dir/"epiJSON_test_1"
python ./run/epiJSON.py \
    --container-image-name "$1" \
    --docker-no-pulling \
    $py_platform_str \
    "--num-threads" "1" \
    "--output-directory" "$epijson_test_dir_1" \
    "--input-format" "JSON" \
    "--input-file" "$dummy_dataset.json" \
    "--phenotype" "DICHOTOMOUS" \
    "--make-all-formats"

cmp $epijson_test_dir_1/dataset.json ${dummy_dataset}_recreated.json
cmp $epijson_test_dir_1/dataset.bim $dummy_dataset.bim
cmp $epijson_test_dir_1/dataset.bed $dummy_dataset.bed
cmp $epijson_test_dir_1/dataset.fam $dummy_dataset.fam
cmp $epijson_test_dir_1/dataset.linden.cases $dummy_dataset.linden.cases
cmp $epijson_test_dir_1/dataset.linden.controls $dummy_dataset.linden.controls
cmp $epijson_test_dir_1/dataset.linden.loci $dummy_dataset.linden.loci
cmp $epijson_test_dir_1/dataset.macoed $dummy_dataset.macoed
cmp $epijson_test_dir_1/dataset.map $dummy_dataset.map
cmp $epijson_test_dir_1/dataset.ped $dummy_dataset.ped
cmp $epijson_test_dir_1/dataset.tfam $dummy_dataset.tfam
cmp $epijson_test_dir_1/dataset.tped $dummy_dataset.tped
cmp $epijson_test_dir_1/dataset.vcf $dummy_dataset.vcf


# TODO: add more tests here
# ...


# remove temp folder
rm -rf ./test_out