#!/usr/bin/env python3
#/////////////////////////////////////////////////////////////////////////////#
#                                                                             #
#   Copyright (C) 2020 by David B. Blumenthal                                 #
#                                                                             #
#   This file is part of NeEDL.                                        #
#                                                                             #
#   NeEDL is free software: you can redistribute it and/or modify it   #
#   under the terms of the GNU General Public License as published by         #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   NeEDL is distributed in the hope that it will be useful,           #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with NeEDL. If not, see <http://www.gnu.org/licenses/>.      #
#                                                                             #
#/////////////////////////////////////////////////////////////////////////////#

##
# @file install.py
# @brief Installs NeEDL and its dependencies.
#
# @details 
# Usage: 
# ```sh
# $ python install.py [--help] [-h] [--target docs|unit|instance|variance_model|penetrance_model|regression_model|bayesian_model|compare_models] [--clean] [--debug]
# ```
#
# For more information, execute `$ python install.py --help`.

'''Installs NeEDL and its dependencies.'''

import subprocess
import argparse
import os.path
import platform

def build_targets(args, all_targets, dependencies):

    if args.clean:
        print("-- Cleaning build directory.")
        commands = "rm -rf build"
        subprocess.call(commands, shell=True)

    print("-- Creating build directory.")
    commands = "mkdir -p build"
    subprocess.call(commands, shell=True)

    if (not os.path.isfile("build/Makefile")):
        print("-- Running CMake.")
        commands = "cd build; rm -rf *; cmake .. -DCMAKE_BUILD_TYPE=" + ("Debug" if args.debug else "Release")

        if args.system_boost:
            commands += f' -DUSE_INCLUDED_BOOST=0'

        commands += f' -DPython_EXECUTABLE="{dependencies["python3"]}"'
        commands += f' -DCMAKE_C_COMPILER="{dependencies["gcc"]}" -DCMAKE_CXX_COMPILER="{dependencies["g++"]}"'
        if platform.system() == "Darwin":
            commands = commands + " -DOMP_ROOT=" + subprocess.check_output("brew --prefix", shell=True).decode("utf-8")
        subprocess.call(commands, shell=True)

    if args.all_targets:
        for target in targets:
            print("-- Building target {}.".format(target))
            commands = "cd build; make {}".format(target)
            subprocess.call(commands, shell=True)
    elif args.target:
        print("-- Building target {}.".format(args.target))
        commands = "cd build; make {}".format(args.target)
        subprocess.call(commands, shell=True)

    print("-- Creating output directories.")
    commands = "mkdir -p test/unit/init; mkdir -p test/model/res"
    subprocess.call(commands, shell=True)

def build_external_libraries(args, dependencies):
    # boost
    if args.system_boost:
        print("-- Boost libraries: system installation is used --> build of boost is skipped")
    elif os.path.isfile("ext/boost_1_71_0/.INSTALLED") and not args.clean:
        print("-- Boost libraries already built.")
    else:
        print("-- Building Boost libraries.")
        # commands = 'cd ext/boost_1_71_0; rm -rf project-config.jam*; ./bootstrap.sh; echo -e "\\nusing mpi ;\\n" >> project-config.jam; ./b2'
        commands = f'cd ext/boost_1_71_0;CC="{dependencies["gcc"]}" CXX="{dependencies["g++"]}";./bootstrap.sh --with-python="{dependencies["python3"]}"'
        subprocess.call(commands, shell=True)
        with open("ext/boost_1_71_0/user-config.jam", 'w') as file:
            file.write(f'using gcc : : {dependencies["g++"]}')
        commands = f'cd ext/boost_1_71_0;CC="{dependencies["gcc"]}" CXX="{dependencies["g++"]}";./b2 --toolset=gcc'
        subprocess.call(commands, shell=True)
        f = open("ext/boost_1_71_0/.INSTALLED", "w")
        f.close()


    # igraph
    if os.path.isfile("ext/igraph_0.9.8/.INSTALLED") and not args.clean:
        print("-- igraph libraries already built.")
    else:
        print("-- Building igraph libraries.")
        commands = f'cd ext/igraph_0.9.8; rm -rf build; mkdir -p build; cd build; cmake .. -DIGRAPH_OPENMP_SUPPORT=on -DIGRAPH_ENABLE_LTO=on -DIGRAPH_ENABLE_TLS=on -DCMAKE_C_COMPILER="{dependencies["gcc"]}" -DCMAKE_CXX_COMPILER="{dependencies["g++"]}"; cmake --build .'
        subprocess.call(commands, shell=True)
        f = open("ext/igraph_0.9.8/.INSTALLED", "w")
        f.close()

    # uWebSockets
    if os.path.isfile("ext/.uWebSockets_INSTALLED") and not args.clean:
        print("-- uWebSocket library already built.")
    else:
        print("-- Building uWebSockets library")
        if not args.no_submodule_extraction:
            subprocess.call('git submodule update --init --recursive')
        
        subprocess.call('cd ext/uWebSockets/uSockets;make boringssl;cd boringssl;BORINGSSL=$PWD;cd ../lsquic;cmake -DBORINGSSL_DIR=$BORINGSSL .;make;cd ..;WITH_LTO=1 WITH_QUIC=1 WITH_BORINGSSL=1 make', shell=True)
        f = open("ext/.uWebSockets_INSTALLED", "w")
        f.close()

def extract_all_zips_in_folder(folder):
    print("extracting files in " + folder)

    if not folder.endswith("/"):
        folder += "/"

    # check if directory even exist --> skip otherwise with warning
    if not os.path.exists(folder):
        print(f"Warning: Path {folder} does not exist. Extracting skipped.")
        return

    # concat all splitted zips
    with os.scandir(folder) as files:
        escaped_folder = folder.replace('"', '\\"')
        for entry in files:
            if entry.name.endswith(".aa"):
                original_name = entry.name[:-3].replace('"', '\\"').replace(" ", '\\ ')
                subprocess.call('cat ' + escaped_folder.replace(" ", '\\ ') + original_name + '.* > ' + escaped_folder.replace(" ", '\\ ') + original_name, shell=True)

    # extract all zip files
    with os.scandir(folder) as files:
        escaped_folder = folder.replace('"', '\\"')
        for entry in files:
            if entry.name.endswith(".zip"):
                print("Extracting " + entry.name)
                subprocess.call('unzip -qo -d "' + escaped_folder + '" "' + escaped_folder + entry.name.replace('"', '\\"') + '"', shell=True)

def extract_resources(args):
    if (os.path.isfile(".RES_EXTRACTED") and not args.clean) or args.no_data_unpacking:
        print("-- data is already extracted.")
    else:
        extract_all_zips_in_folder("ext")
        extract_all_zips_in_folder("data/BIOGRID")
        extract_all_zips_in_folder("data/dbSNP/inc_pseudogenes")
        extract_all_zips_in_folder("data/dbSNP/no_pseudogenes")
        extract_all_zips_in_folder("data/GRN")
        extract_all_zips_in_folder("data/PositionBasedInputModule/snp_balanced_tree")
        extract_all_zips_in_folder("data/PositionBasedInputModule/snp_hashmaps")
        extract_all_zips_in_folder("data/PositionBasedInputModule/gene_balanced_tree")
        extract_all_zips_in_folder("ext/plink")
        f = open(".RES_EXTRACTED", "w")
        f.close()

    if args.extract_datasets:
        if (os.path.isfile(".DS_EXTRACTED") and not args.clean) or args.no_data_unpacking:
            print("-- datasets are already extracted.")
        else:
            extract_all_zips_in_folder("data")


def find_dependencies(args):
    # python executable
    if args.python3 is not None:
        python3 = args.python3
    else:
        commands = "which python3"
        python3 = subprocess.check_output(commands, shell=True, text=True).strip()

    # C/C++ compiler
    if args.gcc is not None:
        gcc = args.gcc
    else:
        commands = "which gcc"
        gcc = subprocess.check_output(commands, shell=True, text=True).strip()

    if args.gxx is not None:
        gxx = args.gxx
    else:
        commands = "which g++"
        gxx = subprocess.check_output(commands, shell=True, text=True).strip()

    print("=====  CONFIGURATION  =====")
    print(f"Python3 Executable: {python3}")
    print(f"gcc: {gcc}")
    print(f"g++: {gxx}")
    print("===========================")

    return {
        "python3": python3,
        "gcc": gcc,
        "g++": gxx
    }


print("\n**************************************************")
print("                       NeEDL                      ")
print("                Installation Script               ")
print("**************************************************")

parser = argparse.ArgumentParser(description="Compiles NeEDL and its dependencies.", epilog="If called without arguments, only the dependencies are installed.")
targets = ["NeEDL", "epiJSON",  "create_JSON_PLINK", "create_JSON_PLINK_bigfiles", "create_LD_matrix", "calculate_scores", "convert_to_binary", "docs", "unit", "instance", "variance_model", "penetrance_model", "regression_model", "bayesian_model", "compare_models", "realtime_scores"]
parser.add_argument("--all-targets", help="build all NeEDL targets", action='store_true')
parser.add_argument("--target", help="build selected target", metavar="NeEDL|epiJSON|create_JSON_PLINK|create_JSON_PLINK_bigfiles|create_LD_matrix|calculate_scores|convert_to_binary|docs|unit|instance|variance_model|penetrance_model|regression_model|bayesian_model|compare_models|realtime_scores", choices=targets)
parser.add_argument("--debug", help="build in debug mode", action="store_true")
parser.add_argument("--clean", help="clean build directory, external libs and update makefile before build", action="store_true")
parser.add_argument("--no-data-unpacking", help="Disables unpacking of data. This option should only be used if the only planned build target is realtime_score.", action="store_true")
parser.add_argument("--python3", help="path to the python3 executable that should be used for building and linking NeEDL", default=None)
parser.add_argument("--gcc", help="one can select the path to gcc manually if the automatically selected compiler is not correct.", default=None)
parser.add_argument("--gxx", help="one can select the path to g++ manually if the automatically selected compiler is not correct.", default=None)
parser.add_argument("--extract-datasets", help="Also extracts simulated datasets. These might be necessary for the unit tests.", action="store_true")
parser.add_argument("--no-submodule-extraction", help="When set the script does not recursively download submodules. This is used for building the docker container.", action="store_true")
parser.add_argument("--system-boost", help="Use the system boost installation instead of the one in this repository", action="store_true")
args = parser.parse_args()


dependencies = find_dependencies(args)

extract_resources(args)

build_external_libraries(args, dependencies)

build_targets(args, targets, dependencies)

print("**************************************************\n")
