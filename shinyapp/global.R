###################
#### Libraries ####
###################

print("######################### LOADING LIBRARIES #########################")

# shiny related packages 
library(shiny)
library(shinyjs)
library(shinydashboard)
library(shinydashboardPlus) 
library(shinyWidgets)
library(shinycssloaders)
library(shinyFeedback)
library(shinyalert)

# customize the appearance and functionality of the R Shiny app 
library(bslib) # for bootstrap
library(rclipboard) # for copy buttons in tables 
library(sortable) # for priority setting of centralities 
library(htmltools) # for HTML generation
library(httr) # for working with HTTP


# table and data cleaning related packages 
library(data.table)
library(DT)
library(tidyr)
library(dplyr)

# for facilitated data handling 
library(stringr)
library(magrittr)


# json handling related packages 
library(rjson)
library(tidyjson) 
library(jsonlite) 


# database releated packages 
library(RSQLite) # for interpreting the sqlite format for the networks


# plotting related packages
library(ggplot2)
library(gghighlight) 
library(plotly)
library(patchwork)

# graph related packages 
library(visNetwork)
library(igraph)


# package for documentation
library(roxygen2)

# package for argument parser
library(argparse)

print("######################### PARSING ARGUMENTS #########################")


parser <- argparse::ArgumentParser()
parser$add_argument("--results", help = "Results directory", required = TRUE)
parser$add_argument("--dataset", help = "Dataset file", required = TRUE)
parser$add_argument("--dbSNP", help = "dbSNP file", required = TRUE)
parser$add_argument("--biogrid", help = "BIOGRID file", required = TRUE)
args <- parser$parse_args()


##########################################
##### REQUIRED INPUT PATHS FROM NEEDL ####
##########################################
# Path to dbSNP 
dbSNP_path <- args$dbSNP

# Path to BIOGRID
biogrid_path <- args$biogrid

# output path with results
out_path <- args$results

print("######################### ARGUMENTS PARSED #########################")


#################
#### Sources ####
#################
print("######################### GET SOURCE FILES #########################")
print("1. TASK STATS")
source("src/create_task_stat.R")

print("2. MAP SNPS TO GENE")
source("src/map_snps_to_genes.R")

print("3. NETWORK")
source("src/network.R")

print("4. CENTRALITIES")
source("src/plot_centralities.R")

print("5. SCORES")
source("src/plot_scores.R")

print("6. LOAD ALL VARIABLES")
source("src/one_run_final.R")




#################
##### RDATA #####
#################
# save_path <- "/nfs/proj/GenEpiSeeker/sylvie/needl_validation/validation_report/r_studio/data/RA/RA_run_biological.RData"
# load(save_path)
# save.image(file = save_path)



######################
#### Session Info ####
######################

# Session info from 2023/05/04

# R version 4.2.1 (2022-06-23)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Ubuntu 22.04.1 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.10.0
# LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.10.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] Rcpp_1.0.10           rstudioapi_0.13       magrittr_2.0.3        bit_4.0.5             xtable_1.8-4         
# [6] rjson_0.2.21          R6_2.5.1              rlang_1.0.6           fastmap_1.1.1         blob_1.2.1           
# [11] tools_4.2.1           shinyWidgets_0.7.6    DT_0.27               DBI_1.1.1             cli_3.5.0            
# [16] jquerylib_0.1.4       htmltools_0.5.4       shinyjs_2.1.0         ellipsis_0.3.2        bit64_4.0.5          
# [21] digest_0.6.31         fontawesome_0.5.0     lifecycle_1.0.3       shiny_1.7.4           later_1.3.0          
# [26] vctrs_0.5.2           htmlwidgets_1.6.2     sass_0.4.4            shinyFeedback_0.4.0   promises_1.2.0.1     
# [31] memoise_2.0.1         RSQLite_2.2.19        cachem_1.0.7          mime_0.12             compiler_4.2.1       
# [36] bslib_0.4.2           shinycssloaders_1.0.0 argparse_2.2.2        jsonlite_1.8.4        httpuv_1.6.7    