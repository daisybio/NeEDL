#!/usr/bin/env python3

import sys
import os

# settings
image = "bigdatainbiomedicine/needl:latest"
command = '/NeEDL/test/model/bin/calculate_scores'

# app arguments
arguments = sys.argv[1:]

print(arguments)


# check if singularity should be used instead of docker
use_singularity = False
singularity_cmd = 'singularity'
for i, arg in enumerate(arguments):
    if arg == '--singularity':
        use_singularity = True
        # remove argument from list as app cannot process it
        del arguments[i]
    if arg == '--apptainer':
        use_singularity = True
        singularity_cmd = 'apptainer'
        # remove argument from list as app cannot process it
        del arguments[i]


# define a different docker image name, needed for testing
for i, arg in enumerate(arguments):
    if arg == '--docker-image-name':
        if len(arguments) > i + 1:
            image = arguments[i + 1]

            # remove arguments from list as NeEDL cannot parse them
            del arguments[i + 1]
            del arguments[i]


# disable pulling of docker container
docker_pull = not use_singularity
for i, arg in enumerate(arguments):
    if arg == '--docker-no-pulling':
        docker_pull = False
        del arguments[i]


# map output path
output_file = None
output_directory = None
for i, arg in enumerate(arguments):
    if arg == '--snp-sets-output-file':
        if len(arguments) > i + 1:
            output_abs = os.path.abspath(arguments[i + 1])
            output_directory = os.path.dirname(output_abs)
            output_file = os.path.basename(output_abs)
            if not os.path.exists(output_directory):
                print(f'Error: output directory {output_directory} does not exist!')
                exit(-1)

            arguments[i + 1] = '/mnt/out/' + output_file

        break


# map all input files and directories
input_paths = []

# normal paths
normal_path_attribs = ['--gwas-input-file', '--snp-sets-input-file']
for i, arg in enumerate(arguments):
    if arg in normal_path_attribs:
        if len(arguments) > i + 1:
            path = os.path.abspath(arguments[i + 1])
            if not os.path.exists(path):
                print(f'Error: path {path} at argument {arg} does not exist!')
                exit(-1)

            arguments[i + 1] = '/mnt/in_' + str(len(input_paths))
            input_paths.append(path)

# special paths
# --snp-annotate ...
for i, arg in enumerate(arguments):
    if arg == '--snp-annotate':
        if len(arguments) > i + 1:
            splits = arguments[i + 1].split('|', maxsplit = 1)
            path = os.path.abspath(splits[0])
            if not os.path.exists(path):
                print(f'Error: path "{path}" at argument "{arg}" does not exist!')
                exit(-1)

            arguments[i + 1] = '/mnt/in_' + str(len(input_paths))
            if len(splits) == 2:
                arguments[i + 1] += '|' + splits[1]

            input_paths.append(path)

# --data-directory
if '--data-directory' in arguments:
    print("Error: It is not allowed to set the argument --data-directory when using NeEDL in docker mode. All necessary files within the ext directory are already at the correct place in the docker container.")
    exit(1)

if use_singularity:
    volume_string = ' '.join([f'-B "{file}:/mnt/in_{i}:rw,Z"' for i, file in enumerate(input_paths)])
    external_command = f"{singularity_cmd} exec {volume_string}"

    if output_directory is not None:
        external_command += f' -B "{output_directory}:/mnt/out:rw,Z" '

    external_command += "docker://"

else:
    volume_string = ' '.join([f'-v "{file}:/mnt/in_{i}:rw,Z"' for i, file in enumerate(input_paths)])
    external_command = f"docker run --user='{os.getuid()}':'{os.getgid()}' {volume_string}"

    if output_directory is not None:
        external_command += f' -v "{output_directory}:/mnt/out:rw,Z" '


argument_string = ' '.join(map(lambda a: f'"{a}"', map(lambda b: b.replace('"', '\\"'), arguments)))
external_command += f"{image} {command} --data-directory /NeEDL/data/ {argument_string}"

print(external_command)

if docker_pull:
    ecode = os.system("docker image pull " + image)
    if ecode != 0:
        sys.exit(1)

ecode = os.system(external_command)
if ecode != 0:
    sys.exit(1)