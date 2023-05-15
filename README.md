# Tools provided in this repository
- NeEDL: Tool for fast network-based epistasis detection using local search
- epiJSON: Tool for converting various input and output formats which also can apply commonly used filtering steps. Needed to create the input files for NeEDL
- calculate_scores: A tool that can be used to calculate score of all statistical models supported by NeEDL (see Blumenthal et al. for further details)

# License and Citing

The project is distributed under the [GNU General Public License](https://www.gnu.org/licenses/gpl-3.0.en.html).

# Usage (Docker container)

**The easiest way to use our toolset is via our docker container**. This repo is big and contains much more than a typical user likely needs. Thus, we strongly encourage you to try out the docker container first, which comes with up-to-date versions of all of our tools in a smaller size.

In order to use our docker container, it is required to have docker installed. The user also needs permission to create and run docker containers (p.r.n., ask your system administrator). Additionally, a Python 3 installation and the command-line tool `curl` are necessary.


```bash
# run NeEDL
curl https://raw.githubusercontent.com/biomedbigdata/NeEDL/main/run/NeEDL.py | python3 - <parameters...>
```
```bash
# run epiJSON
curl https://raw.githubusercontent.com/biomedbigdata/NeEDL/main/run/epiJSON.py | python3 - <parameters...>
```
```bash
# run calculate_scores
curl https://raw.githubusercontent.com/biomedbigdata/NeEDL/main/run/calculate_scores.py | python3 - <parameters...>
```

*todo in einzelne Boxen*

Parameters that can be used with the tools are explained below. Please note the "` - `" before the first parameter. It is necessary to pipe the code of the Python script correctly.

Example command to check that everything works:
```bash
curl https://raw.githubusercontent.com/biomedbigdata/NeEDL/main/run/NeEDL.py | python3 - --help
```


*todo: update mechanism for the container*


## Parameters for NeEDL
### Quick Start / Example
To run NeEDL the default way that we also used during our publication use the following arguments:
```bash
--num-threads <num-threads>
--output-directory <output-directory>
--input-format JSON_EPIGEN
--input-path <input-JSON-file-created-with-epiJSON>
--phenotype DICHOTOMOUS
--snp-annotate-dbSNP
--network-BIOGRID
--ms-seeding-routine RANDOM_CONNECTED
--ms-rc-start-seeds 5000
--ms-model PENETRANCE_NLL
```

The full command then looks like this:
```bash
curl https://raw.githubusercontent.com/biomedbigdata/NeEDL/main/run/NeEDL.py | python3 - \
    --num-threads <num-threads> \
    --output-directory <output-directory> \
    --input-format JSON_EPIGEN \
    --input-path <input-JSON-file-created-with-epiJSON> \
    --phenotype DICHOTOMOUS \
    --snp-annotate-dbSNP \
    --network-BIOGRID \
    --ms-seeding-routine RANDOM_CONNECTED \
    --ms-rc-start-seeds 5000 \
    --ms-model PENETRANCE_NLL
```

Please replace `<...>` with your value of choice.

### Basic Setup
```bash
# if set to 0, all available threads will be used (please set explicit number of threads for use with slurm!)
--num-threads <num-threads>

# Give a path to an existing output directory. Inside the directory, each run will create its unique subdirectory with timestamp and a random number, so multiple runs can use the same output directory. All intermediate and final results (and debug info) are written to this directory.
--output-directory <path>

# Per default, the generated SNP-SNP-interaction network will be stored in the output folder in form of a sqlite database. This database file can get quite big for larger datasets.
# If you do not want this to happen due to memory reasons, please use this flag to deactivate saving the network to disk.
--disable-save-network

# Activates a joint-degree analysis. It creates an additional file containing information about the node degrees of connected SNPs in the SNP-SNP-interaction network. This flag is usually not necessary if one does not want to further investigate the network created by NeEDL.
--do-joint-degree-analysis
```

### Select your dataset
```bash
# Our default format is JSON_EPIGEN. You can use epiJSON to convert various commonly used file formats into this JSON format.
--input-format JSON_EPIGEN

# The actual input file path. Either specify a path to a .json file or path to a directory containing only one .json file.
--input-path <path>

# Phenotype information of the dataset. NeEDL can deal with DICHOTOMOUS, CATEGORICAL and QUANTITATIVE phenotypes. --num-categories is only necessary for CATEGORICAL phenotypes with more than two categories.
--phenotype CATEGORICAL
--num-categories 2
```

### Creation of the SSI network
The creation of the SSI network consists of two phases. First, all SNPs need to be annotated with labels (usually gene labels, but can, in general, be anything). Second, NeEDL creates all edges of the network. SNPs with at least one common label are connected. SNPs that have two different labels are connected if the two labels are set into relation in a label-label interaction network. The default is to use gene names as labels and the protein-protein interaction network BioGRID for the creation of the network.

For this default case which most users probably want to use one has to set the following two flags:
```bash
# annotate SNPs with gene labels from dbSNP
--snp-annotate-dbSNP

# connect SNPs that are labeled with genes which interact according to BioGRID
--network-BIOGRID
```

NeEDL is very versatile and can be used with custom data as well. You can select your data for labeling and interactions:

*tba*

Multiple commands for SNP annotation can be set in a single run to use labels from different sources.

NeEDL also supports using multiple network files in one run. If so, the local search is performed for each network individually and the results are aggregated afterward with the help of an additional local search run. The multiple networks approach is not fully evaluated yet.

### Selecting a seeding routine
NeEDL supports different seeding routines for the local search which come with different advantages and disadvantages. In our publication, we used RANDOM_CONNECTED seeding for all our tests. It is generally a good choice. However, the user has to decide how many seeds should be generated which influences the quality of the results. If you do not know how many seeds to use we suggest trying out the COMMUNITY_WISE seeding as it decides how many seeds should be used based on the network size and architecture.

#### RANDOM_CONNECTED seeding
Random connected seeding selects pairs of SNPs that are connected in the SSI network. In our publication, we used this seeding method with 5000 start seeds for our runs.

```bash
# set the seeding routine
--ms-seeding-routine RANDOM_CONNECTED

# select the number of start seeds to use for the local search. This value has a nearly linear impact on the runtime of NeEDL.
--ms-rc-start-seeds 5000
```
#### COMMUNITY_WISE seeding
Performs community detection on the SSI network producing clusters/communities with a given maximum size. For each cluster, several start seed candidates are selected randomly. Afterward, the actual start seeds are picked from these candidates. The algorithm then makes sure that at least one seed from every cluster is used and that the top $n\%$ of all seed candidates w.r.t. their statistical score is selected.

```bash
# set the seeding routine
--ms-seeding-routine COMMUNITY_WISE

# configure the described options
--ms-cw-quantile 0.25          # pick the 25% best seed candidates
--ms-cw-max-cluster-size 1000  # maximum cluster size allowed
--ms-cw-num-sets-per-cluster 5 # how many seed candidates are generated per cluster
```

#### QUANTUM_COMPUTING seeding
This seeding routine derives from the COMMUNITY_WISE seeing. Hence, all parameters starting with `--ms-cw-...` that are explained above are also available here. Instead of picking random start seeds for each community like the COMMUNITY_WISE seeding does, QUANTUM_COMPUTING seeding uses a quantum computing-based algorithm do select the seed candidates. As this algorithm is potentially resource intensive, one can select to only perform it on clusters of a certain size.

```bash
# set the seeding routine
--ms-seeding-routine QUANTUM_COMPUTING

# hyperparameters

# QC is not applied on every found community but just on communites of a certain minimum size. Smaller communities are handled like COMMUNITY_WISE seeding was selected
--ms-qc-min-cluster-size 100


# Hardware-agnostic hyperparameters
--ms-qc-n-clique 2
--ms-qc-k 3
--ms-qc-nu 3.0
--ms-qc-lambda0 5.0
--ms-qc-lambda1 1.0
--ms-qc-lambda2 1.0
```

The QUANTUM_COMPUTING seeding support three different hardware types. Options are SIMULATED_ANNEALING, QUANTUM_ANNEALING, and QAOA. Depending on the selected hardware different additional parameters are needed as you can see below:

```bash
# for SIMULATED_ANNEALING
--ms-qc-mode SIMULATED_ANNEALING
--ms-qc-sa-num-samples 10
--ms-qc-sa-num-sweeps 1000
--ms-qc-sa-seed 1234


# for QUANTUM_ANNEALING
--ms-qc-mode QUANTUM_ANNEALING
--ms-qc-qa-token "abc"
--ms-qc-qa-num-reads 10
--ms-qc-qa-solver-idx 4
--ms-qc-qa-fw-annealing-ramp-time 1.0
--ms-qc-qa-fw-annealing-pause-time 1.0
--ms-qc-qa-rev-annealing-ramp-time 1.0
--ms-qc-qa-rev-annealing-pause-time 1.0
--ms-qc-qa-rev-annealing-s-target 1.0

# for QAOA
--ms-qc-mode QAOA
--ms-qc-qaoa-vendor abc
--ms-qc-qaoa-azure-subscription-id abc
--ms-qc-qaoa-azure-location abc
--ms-qc-qaoa-azure-backend abc
--ms-qc-qaoa-optimizer abc
--ms-qc-qaoa-maxiter 10
--ms-qc-qaoa-reps 5
--ms-qc-qaoa-n-shots 5
--ms-qc-qaoa-is-recursive-qaoa 0
```

### Local search
One can influence how the local search that NeEDL applies to the SSI network behaves. The most important setting is the statistical score that should be used for optimization. We recommend to use the score ```PENETRANCE_NLL``` for optimization.
```bash
# select score for optimization
--ms-model PENETRANCE_NLL

# optionally select a time limit for search
--ms-per-seed-time-limit 10m   # default is no limit
--ms-per-search-time-limit 3d  # default is no limit

# The arguments below usually are not necessary to be set. The change the behaviour of the local search and simulated annealing.
--ms-max-rounds 300
--ms-min-set 2 
--ms-max-set 10    
--ms-annealing-type SIMULATED_ANNEALING
--ms-cooling-factor 1.0
--ms-annealing-start-prob 0.8
--ms-annealing-end-prob 0.01

```

*todo: explain additional parameters*

## Documentation of the output files of NeEDl

*tba*

## Parameters for epiJSON
epiJSON is a wrapper around plink that can convert several plink-readable input files and apply several filters typically used for GWAS. It can also create the JSON_EPIGEN files that NeEDL requires as input. The following input and output formats are supported:
- `.bim/.bed./.fam`
- `.ped/.map`
- `.tped/.tfam`
- `.vcf`
- `.bcf2`
- `.json`: JSON_EPIGEN (input for NeEDL)
- MACOED input files (*can only be generated not read by epiJSON*)
- LINDEN input files (*can only be generated not read by epiJSON*)

### Quick Start
Quickly convert your files into the correct format for NeEDL without applying any optional filters:
```bash
curl https://raw.githubusercontent.com/biomedbigdata/NeEDL/main/run/epiJSON.py | python3 - \
    --num-threads <num-threads> \
    --output-directory <output-directory> \
    --input-file <path-to-input> \
    --input-format <BIM_BED_FAM|PED_MAP|TPED_TFAM|VCF|BCF2> \
    --phenotype <DICHOTOMOUS|CATEGORICAL|QUANTIATIVE> \
    --make-json
```

### General parameters
```bash
# Specifiy how many threads epiJSON should use. Most steps will still be run single-threaded.
--num-threads <num-threads>

# Specify an empty and existing output directory where all final files will be written to
--output-directory <path-to-existing-empty-directory>

# Set the phenotype
--phenotype <DICHOTOMOUS|CATEGORICAL|QUANTITATIVE>

# Optionally set the number of categories for a CATEGORICAL phenotype
--num-categories <number>
```

### Selecting the input
epiJSON supports joining multiple input files (possibly of different file type) together. This can be especially helpful if one has cases and controls in separate files. epiJSON can also override the phenotype for each input file.

For every input file specify these parameters:
```bash
--input-file <path>

# Explicitly tell epiJSON the input format. If omitted, epiJSON tries to guess which file to use.
--input-format <BIM_BED_FAM|PED_MAP|TPED_TFAM|VCF|BCF2>

# By default epiJSON uses the phenotype reported in the input file. Optionally, one can specify a phenotype for each file that should be used instead of the reported phenotype that is contained in the input.
--override-phenotype <value>
```

Example for combining separate case and control files:
```bash
curl https://raw.githubusercontent.com/biomedbigdata/NeEDL/main/run/epiJSON.py | python3 - \
    --num-threads 1 \
    --output-directory combined_dataset \
    \
    --input-file cases \
    --input-format BIM_BED_FAM \
    --override-phenotype 2 \
    \
    --input-file controls \
    --input-format BIM_BED_FAM \
    --override-phenotype 1 \
    \
    --phenotype DICHOTOMOUS \
    --make-json
```

`--override-phenotype` and `--input-format` can be specified only once if the same value should be used for all input files.

There is the special case where one might want to override the phenotype only for some input files but not all. To do that, please set `--override-phenotype` to `NO` for the files that you do not want to override. If all phenotype values should be sourced from the input file, then simply do not specify any `--override-phenotype` parameter.

### Filter selection
```bash
# select all available filters
--all-filters


# Filter out any variants which IDs do not start with 'rs'
--filter-rsids

# Apply plink's duplicate filter
--filter-duplicates

# Filters out all genotypes and individuals with > 0.2 missing values.
--filter-missing

# Filters out samples with discrepancy between reported sex and sex according to genome.
--filter-sex-discrepance

# Apply MAF filter with threshold 0.05 for sets with #individuals < 50,000 and 0.01 else.
--filter-allele-frequency

# Filter based on the Hard-Weinberg test.
--filter-hardy-weinberg

# Erases heterzygous haploid and nonnormal Y chromosome genotype calls.
--filter-hh-missing

# Filters out all SNPs that would not be part of the dbSNP/BioGRID based SNP-SNP-interaction network used by NeEDL.
--filter-snp-network
```

### Output files
```bash
# Generate all available output files. (LINDEN and MACOED files will only be created if phenotype is DICHOTOMOUS)
--make-all-formats

# Create JSON_EPIGEN for NeEDL
--make-json

# Create .bim/.bed/.fam files
--make-bim-bed-fam

# Create .ped/.map files
--make-ped-map

# Create .tped/.tfam files
--make-tped-tfam

# Create a VCFv4.2 file which is block-gzipped
--make-vcf

# Create input files for the tool LINDEN
--make-linden

# Create input file for the tool MACOED
--make-macoed
```

## Parameters for calculate_scores

*todo*


# Usage (via repository)
We do not recommend building the tools in this project yourself as our docker container is easier and does not require the user to install any additional dependencies. The following documentation is primarily meant for users who want to extend or modify our code and thus are required to build from scratch. If you are just a normal user and the docker container does not work on your platform please do not hesitate to create a github issue so that we can extend compatibilty.

## Dependencies

Before using GenEpiSeeker, you have to install the following external dependencies:

- [CMake](https://cmake.org/) version 2.6 or higher (for compilation). Installation instructions can be found [here](https://cmake.org/install/). 
- [Doxygen](http://www.stack.nl/~dimitri/doxygen/) (for creating the documentation). Installation instructions can be found [here](https://www.stack.nl/~dimitri/doxygen/manual/install.html).
- [OpenMP](http://www.openmp.org/) compatible C++ compiler. Under Linux, OpenMP is supported by default. Under macOS, please install [libomp](https://formulae.brew.sh/formula/libomp) using [Homebrew](https://brew.sh/).  After having installed Homebrew, open a shell and execute `brew install libomp`.

The following external dependencies are distributed with GenEpiSeeker and do not have to be installed separately:

- [Boost (version 1.71.0)](https://www.boost.org/).
- [Catch (version 2.13.9)](https://github.com/catchorg/Catch2)
- [CLI11 (version 1.9.0)](https://github.com/CLIUtils/CLI11)
- [Eigen (version 3.3.7)](http://eigen.tuxfamily.org/)

*TODO: check dependencies*

## Building the tools

The first time, `install.py` will also compile and set up the dependencies in the `ext/` and `data/` folders.

```bash
# clone this repo
git clone https://github.com/biomedbigdata/NeEDL.git NeEDL
cd NeEDL

# get usage and available build targets
./install.py --help

# build only NeEDL (this sould be sufficient for most users)
./install.py --target NeEDL

# alternatively build all targets (might need additional dependencies)
./install.py --all-targets
# or build a different target
./install.py --target <target>
```

If the installation routine does not find the (correct) location for Python or the gcc compiler, one can specify the paths manually:
```bash
./install.py ... --python3 <path> --gcc <path> --gxx <path>
```

Per default, the dependencies are only compiled when using `./install.py` the first time. To recompile the dependencies before the build use the following flag:
```bash
.install.py ... --clean
```

The final binaries can be found at `./test/model/bin`.

## Additional parameters
All of the parameters documented in the section **Usage (Docker container)** can be used analogously. Some additional parameters are necessary when using one of the tools outside the docker container. The parameters mentioned in this section are not necessary (and not allowed) to be used with the docker container.

```bash
# for NeEDL, calculate_scores and epiJSON:
# Either run NeEDL with cwd set to the root of the repo or pass this parameter specifing the path to the data/ directory of the repo (must end with a /).
--data-directory <path>

# for epiJSON:
# Either run NeEDL with cwd set to the root of the repo or pass this parameter specifing the path to the ext/ directory of the repo (must end with a /).
--ext-directory <path>
```
