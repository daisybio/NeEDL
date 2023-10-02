#!/bin/bash

# This scripts tests the functionality of the docker container
# Please provide the docerk image name as first parameter to this script.


# ---- setup

# abort on first error
set -e

# die method to fail the test
warn () {
    echo "$0:" "$@" >&2
}
die () {
    rc=$1
    shift
    warn "$@"
    exit $rc
}

# create temp output folder
mkdir -p ./test_out
out_dir=$(realpath ./test_out)

echo "Temp output directory: $out_dir"
echo "Working directory: $(pwd)"


# ---- tests

# check if binaries exist and are executable
python ./run/NeEDL.py --docker-image-name "$1" --docker-no-pulling --help
python ./run/epiJSON.py --docker-image-name "$1" --docker-no-pulling  --help
python ./run/calculate_scores.py --docker-image-name "$1" --docker-no-pulling  --help
python ./run/convert_to_binary.py --docker-image-name "$1" --docker-no-pulling  --help

# realtime_scores is a special case --> no python launcher script exists for it
# run containers as current user
user_id=$(id -u)
group_id=$(id -g)
docker run --user $user_id:$group_id "$1" /NeEDL/test/model/bin/realtime_scores --help



# select small dummy dataset
dummy_dataset=$(realpath ./data/e2e_tests/dummy_dataset.json)

# test NeEDL quick start example with dummy dataset
python ./run/NeEDL.py \
    --docker-image-name "$1" \
    --docker-no-pulling \
    --num-threads 1 \
    --output-directory "$out_dir" \
    --input-format JSON_EPIGEN \
    --input-path "$dummy_dataset" \
    --phenotype DICHOTOMOUS \
    --snp-annotate-dbSNP \
    --network-BIOGRID \
    --ms-seeding-routine RANDOM_CONNECTED \
    --ms-rc-start-seeds 5 \
    --ms-model PENETRANCE_NLL


# test that custom annotation files and custom networks work
annotation_file=$(realpath ./data/e2e_tests/annotations.csv)
network_file=$(realpath ./data/e2e_tests/network.csv)

# currently this creates an empty network --> checks only if annotation and network file is processed correctly
# Maybe replace with something that actually creates a network from the custom data?
python ./run/NeEDL.py \
    --docker-image-name "$1" \
    --docker-no-pulling \
    --num-threads 1 \
    --output-directory "$out_dir" \
    --input-format JSON_EPIGEN \
    --input-path "$dummy_dataset" \
    --phenotype DICHOTOMOUS \
    --snp-annotate "$annotation_file|yes|SNP|gene|\t|-1|-1" \
    --network "dummy_network|$network_file|yes|col1|col2|\t|-1|-1" \
    --ms-seeding-routine RANDOM_CONNECTED \
    --ms-rc-start-seeds 5 \
    --ms-model PENETRANCE_NLL
 
covariates_file=$(realpath ./data/COV_TEST/EpiGEN_RND_COV.csv)

python ./run/NeEDL.py \
    --input-path "$dummy_dataset" \
    --input-format JSON_EPIGEN \
    --phenotype DICHOTOMOUS \
    --num-threads 1 \
    --output-directory "$out_dir" \
    --snp-annotate-dbSNP \
    --network-BIOGRID \
    --ms-model PENETRANCE \
    --ms-seeding-routine RANDOM_CONNECTED \
    --ms-rc-start-seeds 5 \
    --ms-search-time-limit 5s \
    --covariates-file "$covariates_file"

# check that no NA values are present
for dir in "$out_dir"/* ; do
    if [ -f "$dir/BIOGRID_seeds.csv" ]; then
        if grep -Fxq "nan" "$dir/BIOGRID_seeds.csv"; then
            die 1 "'nan' found in $dir/BIOGRID_seeds.csv"
        fi
    fi

    if [ -f "$dir/BIOGRID_results.csv" ]; then
        if grep -Fxq "nan" "$dir/BIOGRID_results.csv"; then
            die 1 "'nan' found in $dir/BIOGRID_results.csv"
        fi
    fi

done


# TODO: add more tests here
# ...


# remove temp folder
rm -rf ./test_out