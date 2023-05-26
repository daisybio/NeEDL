##########################
#### One run of NeEDL ####
##########################


########################################
##### Structure for output of NeEDL ####
########################################

# File with score for each SNP in a candidate SNP set: BIOGRID_ind_SNP_scores.csv
# SNP-SNP interaction network: BIOGRID_network.sqlite3
# Candidate SNP sets with rs_ids and scores: BIOGRID_results.csv
# Search score over time: BIOGRID_search_score_over_time.csv
# Start seeds of BIOGRID: BIOGRID_seeds.csv
# Configuration of pipeline: pipeline_config.json
# Log file of NeEDL run: run.log


##################
##### General ####
##################

addResourcePath(prefix = "assets", directoryPath = "www")

# About NeEDL
# Image of NeEDL
needl_img_path <- "assets/NeEDL_static_border.png"

# Workflow of NeEDL
# Image of NeEDL workflow
needl_workflow_img_path <- "assets/NeEDL_complete_overview_incLayers_18diseases.png"


############################
##### Technical Details ####
############################

# Parameter information about NeEDL run from pipeline_config.json
pipeline_config_json <- rjson::fromJSON(file = paste0(out_path, "pipeline_config.json"))
pipeline_config_dt <- create_param_from_pipeline_config(pipeline_config_json)


# Log file of NeEDL run: run.log
run_log_dt <- create_stats_from_run_log(paste0(out_path, "run.log"))


# Runtime (from run log file)

# Memory (from run log file)



#####################
##### Statistics ####
#####################

# Start seeds of BIOGRID: BIOGRID_seeds.csv
seeds_dt <- fread(paste0(out_path, "BIOGRID_seeds.csv"))

# Candidate SNPs with individual scores BIOGRID_ind_SNP_scores.csv
snp_scores_dt <- fread(paste0(out_path, "BIOGRID_ind_SNP_scores.csv"))
snp_scores_dt$RS_IDS <- gsub("rs", "", snp_scores_dt$RS_IDS)


# Candidate SNP sets with rs_ids and scores
res_ext_dt <- fread(paste0(out_path, "BIOGRID_results.csv"))


###############################

# Search score over time
search_score_dt <- fread(paste0(out_path, "BIOGRID_search_score_over_time.csv"))


# Set scores
score_names <- colnames(res_ext_dt)[grepl("PENETRANCE|BAYESIAN|REGRESSION|VARIANCE", colnames(res_ext_dt))]
score_names <- unique(score_names[!grepl("RANK", score_names)])



#############################
#### NETWORK INFORMATION ####
#############################

# SNP-SNP interaction network
net_path <- paste0(out_path, "BIOGRID_network.sqlite3")


### Create ppi edges, nodes, i_nodes, i_edges, i_graph for centrality measures
con <- sqlite_load(net_path)
edges <- sqlite_get_dt(con, "edges")
nodes <- sqlite_get_dt(con, "nodes")
colnames(nodes)[2] <- "rs_id"

# Nodes and edges for igraph
i_edges <- get_i_edges(edges)
i_nodes <- get_i_nodes(nodes)

# igraph object
i_graph <- get_igraph_obj(i_nodes, i_edges)


### Network centralities
nodes <- get_cand_centr(i_graph, nodes, nodes$rs_id, centr_measure = "degree", mode = "all")
nodes <- get_cand_centr(i_graph, nodes, nodes$rs_id, centr_measure = "eigenvector", mode = "all")
nodes <- get_cand_centr(i_graph, nodes, nodes$rs_id, centr_measure = "pagerank", mode = "all")




# Take the SNPs of the first n candidate sets to calculate their closeness
first_n_cands <- c(unique(unlist(strsplit(res_ext_dt[, RS_IDS], ";")))[1:100], nodes[, rs_id])[1:1000]

nodes <- get_cand_centr(i_graph, nodes, first_n_cands, centr_measure = "closeness", mode = "all")



# Remove rs from rsIDs
nodes$rs_id <- gsub("rs", "", nodes$rs_id)

# Add scores to nodes
nodes <- merge(nodes, snp_scores_dt, by.x = "rs_id", by.y = "RS_IDS", all = TRUE)



##################################
##### Epistasis disease atlas ####
##################################

# Image of epistasis disease atlas
epiatlas_img_path <- "assets/Epistasis_static_border.png"

# Information for SNP-gene-mapping
dbSNP_dt <- fread(dbSNP_path, header = FALSE)
colnames(dbSNP_dt) <- c("gene", "SNPs")

# Information about protein-protein interaction
biogrid_dt <- fread(biogrid_path)


## Databases for SNPs

# dbSNP
dbSNP_img_path <- "assets/dbSNP_broad.PNG"

# ClinVar
clinvar_img_path <- "assets/ClinVar_broad.PNG"

# SNPmap
snpmap_img_path <- "assets/SNPmap.PNG"


## Databases for genes

# Genecards
genecards_img_path <- "assets/GeneCards_text.png"

# Genbank
genbank_img_path <- "assets/GenBank_broad.PNG"

# Wikipathways
wikipathways_img_path <- "assets/WikiPathways_text.png"

# Gene ontology
geneontology_img_path <- "assets/GeneOntology_text.png"

# Cytoscape
cytoscape_img_path <- "assets/Cytoscape_text.PNG"


## Gene enrichment with gprofiler
gprofiler_img_path <- "assets/gprofiler.png"
