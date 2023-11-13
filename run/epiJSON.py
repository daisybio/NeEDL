#!/usr/bin/env python3

import sys
import os
from pathlib import Path

# settings
image = "bigdatainbiomedicine/needl:latest"
command = '/NeEDL/test/model/bin/epiJSON'

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


# on selinux systems the Z flag might be necessary to mount the directories correctly to the docker container
docker_selinux_flag = ''
for i, arg in enumerate(arguments):
    if arg == '--docker-selinux':
        docker_selinux_flag = ',Z'
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
output_directory = None
for i, arg in enumerate(arguments):
    if arg == '--output-directory':
        if len(arguments) > i + 1:
            output_directory = os.path.abspath(arguments[i + 1])
            if not os.path.exists(output_directory):
                print(f'Error: output directory {output_directory} does not exist!')
                exit(-1)

            arguments[i + 1] = '/mnt/out/'

        break

# map all input files and directories
input_paths = []

# normal paths
normal_path_attribs = ['--input-file', '--disease-snps']
for i, arg in enumerate(arguments):
    if arg in normal_path_attribs:
        if len(arguments) > i + 1:
            # can be multiple files (plink) so always mount the whole directory
            path = os.path.abspath(arguments[i + 1])
            filename = ""
            if not os.path.isdir(path):
                filename = os.path.basename(path)
                path = str(Path(path).parent.resolve())
            if not os.path.exists(path):
                print(f'Error: path {path} at argument {arg} does not exist!')
                exit(-1)

            arguments[i + 1] = '/mnt/in_' + str(len(input_paths)) + '/' + filename
            input_paths.append(path)


# --ext-directory
if '--ext-directory' in arguments:
    print("Error: It is not allowed to set the argument --ext-directory when using NeEDL in docker mode. All necessary files within the ext directory are already at the correct place in the docker container.")
    exit(1)

# --data-directory
if '--data-directory' in arguments:
    print("Error: It is not allowed to set the argument --data-directory when using NeEDL in docker mode. All necessary files within the data directory are already at the correct place in the docker container.")
    exit(1)

if use_singularity:
    volume_string = ' '.join([f'-B "{file}:/mnt/in_{i}:rw"' for i, file in enumerate(input_paths)])
    external_command = f"{singularity_cmd} exec {volume_string}"

    if output_directory is not None:
        external_command += f' -B "{output_directory}:/mnt/out:rw" '

    external_command += "docker://"

else:
    volume_string = ' '.join([f'-v "{file}:/mnt/in_{i}:rw{docker_selinux_flag}"' for i, file in enumerate(input_paths)])
    external_command = f"docker run --user='{os.getuid()}':'{os.getgid()}' {volume_string}"

    if output_directory is not None:
        external_command += f' -v "{output_directory}:/mnt/out:rw{docker_selinux_flag}" '


argument_string = ' '.join(map(lambda a: f'"{a}"', map(lambda b: b.replace('"', '\\"'), arguments)))
external_command += f"{image} {command} --ext-directory /NeEDL/ext/ --data-directory /NeEDL/data/ {argument_string}"

print(external_command)

if docker_pull:
    ecode = os.system("docker image pull " + image)
    if ecode != 0:
        sys.exit(1)

ecode = os.system(external_command)
if ecode != 0:
    sys.exit(1)