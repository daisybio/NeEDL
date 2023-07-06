#!/bin/bash

# This scripts tests the functionality of the docker container
# Please provide the docerk image name as first parameter to this script.


# ---- setup

# abort on first error
set -e

# run containers as current user
user_id=$(id -u)
group_id=$(id -g)


# ---- tests

# check if binaries exist and are executable
docker run --user $user_id:$group_id "$1" /NeEDL/test/model/bin/NeEDL --help
docker run --user $user_id:$group_id "$1" /NeEDL/test/model/bin/epiJSON --help
docker run --user $user_id:$group_id "$1" /NeEDL/test/model/bin/calculate_scores --help
docker run --user $user_id:$group_id "$1" /NeEDL/test/model/bin/convert_to_binary --help


# TODO: add more tests here
# ...