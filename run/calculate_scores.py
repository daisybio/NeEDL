#!/usr/bin/env python3

import sys
import os

# settings
image = "bigdatainbiomedicine/needl:latest"
command = '/NeEDL/test/model/bin/calculate_scores'

# app arguments
arguments = sys.argv[1:]

print(arguments)

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
            splits = arguments[i + 1].split('|', maxsplit = 2)
            path = os.path.abspath(splits[0])
            if not os.path.exists(path):
                print(f'Error: path "{path}" at argument "{arg}" does not exist!')
                exit(-1)

            arguments[i + 1] = '/mnt/in_' + str(len(input_paths))
            input_paths.append(path)

# --data-directory
if '--data-directory' in arguments:
    print("Error: It is not allowed to set the argument --data-directory when using NeEDL in docker mode. All necessary files within the ext directory are already at the correct place in the docker container.")
    exit(1)

volume_string = ' '.join([f'-v "{file}:/mnt/in_{i}:rw,Z"' for i, file in enumerate(input_paths)])

external_command = f"docker run --user='{os.getuid()}':'{os.getgid()}' {volume_string}"

if output_directory is not None:
    external_command += f' -v "{output_directory}:/mnt/out:rw,Z" '

argument_string = ' '.join(map(lambda a: f'"{a}"', map(lambda b: b.replace('"', '\\"'), arguments)))

external_command += f"{image} {command} --data-directory /NeEDL/data/ {argument_string}"

print(external_command)

# os.system("docker image pull " + image)
os.system(external_command)
