#!/bin/bash

# This scripts tests the functionality of the docker container
# Please provide the docerk image name as first parameter to this script.


# ---- setup

# abort on first error
set -e

# run containers as current user
user_id=$(id -u)
group_id=$(id -g)

# create temp output folder
mkdir -p ./test_out
out_dir=$(realpath ./test_out)

echo "Temp output directory: $out_dir"
echo "Working directory: $(pwd)"


# ---- tests

# check if binaries exist and are executable
docker run --user $user_id:$group_id "$1" /NeEDL/test/model/bin/NeEDL --help
docker run --user $user_id:$group_id "$1" /NeEDL/test/model/bin/epiJSON --help
docker run --user $user_id:$group_id "$1" /NeEDL/test/model/bin/calculate_scores --help
docker run --user $user_id:$group_id "$1" /NeEDL/test/model/bin/convert_to_binary --help


# test NeEDL quick start example with dummy dataset
docker run --user $user_id:$group_id \
    -v "./test/e2e/dummy_dataset.json:/mnt/in_0:rw,Z" \
    -v "$out_dir:/mnt/out:rw,Z" \
    "$1" \
    /NeEDL/test/model/bin/NeEDL \
    --data-directory /NeEDL/data/ \
    "--num-threads" "1" \
    "--output-directory" "/mnt/out" \
    "--input-format" "JSON_EPIGEN" \
    "--input-path" "/mnt/in_0" \
    "--phenotype" "DICHOTOMOUS" \
    "--snp-annotate-dbSNP" \
    "--network-BIOGRID" \
    "--ms-seeding-routine" "RANDOM_CONNECTED" \
    "--ms-rc-start-seeds" "5" \
    "--ms-model" "PENETRANCE_NLL"


# TODO: add more tests here
# ...


# remove temp folder
rm -rf ./test_out