##########################
#### Network analysis ####
##########################


#################
#### sqlite #####
#################
# Function to connect to db with dbname = "... .sqlite3"
sqlite_load <- function(dbname){
  con <- dbConnect(drv=RSQLite::SQLite(), dbname)
}

# Function to get structure of sqlite object
sqlite_str <- function(con){
  tables <- dbListTables(con)
  print("TABLES")
  print(tables)
}

# Function to get structure of table in sqlite object
sqlite_str_table <- function(con, table){
  print("FIELDS")
  print(dbListFields(con, table))
}

# Function to create data tables from the sqlite tables 
sqlite_to_dts <- function(con){
  tables <- dbListTables(con)
  # exclude sqlite_sequence (contains table information)
  tables <- tables[tables != "sqlite_sequence"]
  dt_list <- vector("list", length=length(tables))
  ## create a data.frame for each table
  for (i in seq(along=tables)) {
    dt_list[[i]] <- data.table(dbGetQuery(conn=con, statement=paste("SELECT * FROM '", tables[[i]], "'", sep="")))
  }
  
  return(dt_list)
}

# Function to create data tables from the sqlite tables 
# 1: "edges", 2: "has_annotation", 3: "node_annotations", 4: "nodes" 
sqlite_get_dt <- function(con, table){

    table_dt <- data.table(dbGetQuery(conn=con, statement=paste0("SELECT * FROM '", table, "'")))
  
  return(table_dt)
}


################
#### igraph ####
################

# Function to get nodes from dt list
get_i_nodes <- function(nodes){
  i_nodes <- data.frame(nodes[, .(id)])
  
  return(i_nodes)
}

# Function to get edges from dt list 
get_i_edges <- function(edges){
  i_edges <- data.frame(edges[, .(node1, node2)])
  colnames(i_edges) <- c("from", "to")
  
  return(i_edges)
}

# Function to get igraph object from nodes and edges 
get_igraph_obj <- function(i_nodes, i_edges, directed = FALSE){
  igraph_obj <- graph_from_data_frame(i_edges, directed, i_nodes)
  
  return(igraph_obj)
}


#############################
#### Centrality measures ####
#############################


#' Function to get required centrality for candidate SNPs 
#'
#' @param i_graph I_graph object of the network
#' @param nodes Nodes data table with id and rs_id 
#' @param cands All candidate rs_ids for which the centr_measure should be calculated (default: all nodes)
#' @param centr_measure "degree", "closeness", "betweeness", "eigenvector", "pagerank"
#' @param mode "out” for out-degree, “in” for in-degree or “all” for both (default: all)
#' @param directed Indicating if network has directed/undirected edges (default: undirected)
#' @param cutoff Cutoff for betweenness (default: 3) 
#'
#' @return Data table with id, rs_id, calculated centr_measure 
get_cand_centr <- function(i_graph, nodes, cands = nodes$rs_id, centr_measure = "degree", mode = "all", directed = FALSE, cutoff = 3){
  # Get candidate ids using nodes data table
  cand_ids <- nodes[rs_id %in% cands, id]
  cand_ids <- cand_ids + 1
  
  # Create data table with ids and centrality measure to merge with nodes data table
  # Centrality column
  centr_col = c()
  if (centr_measure == "degree"){
    centr_col <- degree(i_graph, V(i_graph)[cand_ids], mode)
  } else if (centr_measure == "closeness"){
    centr_col <- closeness(i_graph, V(i_graph)[cand_ids], mode)
  } else if (centr_measure == "betweenness"){
    centr_col <- betweenness.estimate(i_graph, V(i_graph)[cand_ids], directed, cutoff = cutoff)
  } else if (centr_measure == "eigenvector"){
    centr_col <- eigen_centrality(i_graph, directed)
    centr_col <- centr_col$vector[cand_ids]
  } else if (centr_measure == "pagerank"){
    centr_col <- page_rank(i_graph, algo = "prpack", V(i_graph)[cand_ids], directed)
    centr_col <- centr_col$vector
  }
  
  
  # Id column  
  id_col <- as.numeric(names(centr_col))
  
  # Id centrality data table 
  id_centr_dt <- data.table(id_col, centr_col)
  colnames(id_centr_dt) <- c("id", eval(centr_measure))

  # Merge id centrality data table with nodes data table 
  merge_dt <- merge(nodes, id_centr_dt, by = c("id"), all = TRUE)
  
  
  return(merge_dt)
}



#' Function to sort a network according to their centralities
#'
#' @param rs_id_centr_dt Data table with columns: id, rs_id, centrality/centralities 
#' @param centr_filter Vector containing strings by which the data table will be filtered (default: degree)
#' @param ascending Vector with -1 (descending) and 1 (ascending) (default: -1)
#' @param num_res Number of results to return (default: 10)
#'
#' @return Data table with the first num_res rows filtered by centr_filter, ascending 

network_filter <- function(rs_id_centr_dt, centr_filter = c("degree"), ascending = c(-1), num_res = 10){
  rs_id_centr_dt_copy <- copy(rs_id_centr_dt)
  setorderv(rs_id_centr_dt_copy, centr_filter, ascending)
  return(rs_id_centr_dt_copy[0:num_res])
}



#############################
#### Network statistics  ####
#############################

#' Function to get the statistics of a network e.g., #nodes, #edges
#'
#' @param nodes Nodes of the network
#' @param edges Edges of the network
#'
#' @return Data table with network statistics 
#' @example   
#' get_net_stats(nodes, edges)
get_net_stats <- function(nodes, edges){
  net_stats_dt <- data.table(
    "#nodes" = dim(nodes)[1],
    "#edges" = dim(edges)[1], 
    "#edges SAME TAG" = sum(edges$SAME_TAG), 
    "#edges BIOGRID" = sum(edges$BIOGRID)
  )
  
  colnames_saved <- colnames(net_stats_dt)
  
  net_stats_dt <- t(net_stats_dt)
  net_stats_dt <- data.table(
    Statistic = colnames_saved,
    Value = net_stats_dt[, c(1)]
    )
  
  return(net_stats_dt)
}




###############################
#### Network Visualization ####
###############################


#' Function to plot a network by taking the desired node ids and edges
#'
#' @param rs_id_centr_dt Data table with columns: id, rs_id, centrality/centralities 
#' @param edges Edges of the network
#' @param rs_ids rs ids of nodes on which the network will be based on 
#' @param filter Boolean value with "TRUE" for filtering, "FALSE" otherwise (default: TRUE)
#' @param centr_filter Vector containing strings by which the data table will be filtered (default: degree)
#' @param ascending Vector with -1 (descending) and 1 (ascending) (default: -1)
#' @param num_res Number of results to return (default: 10)
#'
#' @return List with SNP-SNP interaction network, nodes and edges 
#'
#' @examples
#' rs_ids <- c("rs7596121", "rs41440544", "rs41501252")
#'example_net <- show_network_from_id(rs_id_centr_dt = nodes, edges, rs_ids, filter = TRUE,
#'                                    centr_filter = c("degree"), ascending = c(-1), num_res = 10)
show_network_from_id <- function(rs_id_centr_dt, edges, rs_ids, filter = TRUE, 
                                 centr_filter = c("degree"), ascending = c(-1), num_res = 10){
  
  # Substitute rs in the rs_ids by ""
  node_ids <- gsub("rs", "", rs_ids)
  # Get all ids for the given rs_ids
  node_ids <- rs_id_centr_dt[rs_id %in% node_ids, id]
  
  # nodes only with ids 
  nodes_viz <- data.frame(id = node_ids)
  
  # edges from node ids
  # edges from node ids
  edges_viz <- data.frame(edges[node1 %in% nodes_viz$id | node2 %in% nodes_viz$id])
  colnames(edges_viz)[1] <- "from"
  colnames(edges_viz)[2] <- "to"
  
  # nodes extended with all nodes from edges
  nodes_viz <- data.table(id = unique(c(edges_viz$from, edges_viz$to)))
  
  
  # Filtering by centralities using function network_filter
  if(filter){
    nodes_viz <- merge(rs_id_centr_dt, nodes_viz, by = "id")
    
    # Set num_res in case the number of entries in rs_id_centr_dt is lower
    nodes_viz_wo_cands <- nodes_viz[!(id %in% node_ids)]
    
    num_res <- ifelse(dim(nodes_viz_wo_cands)[1] < num_res, dim(nodes_viz_wo_cands)[1], num_res)
    
    # Use network filter function 
    nodes_viz <- network_filter(nodes_viz_wo_cands, centr_filter, ascending, num_res)
    # Create data table with node_ids and the nodes_viz from the filtering step
    nodes_viz <- data.table(id = c(nodes_viz$id, node_ids))
  }
  
  ### Nodes 
  # Get unique node ids ordered by id
  nodes_viz <- data.frame(unique(nodes_viz)[order(id)])
  
  nodes_viz <- merge(rs_id_centr_dt, nodes_viz, by = "id")
  
  
  # label
  colnames(nodes_viz)[2] <- "label"
  nodes_viz$label <- paste0("rs", nodes_viz$label)
  
  # label color 
  font_color_col <- rep("black", dim(nodes_viz)[1])
  nodes_viz[, font.color := font_color_col] 
  
  
  # color
  color_col <- ifelse(nodes_viz$id %in% node_ids, "orange", "lightskyblue")
  nodes_viz[, color.background := color_col] 
  
  # border color
  border_color_col <- ifelse(nodes_viz$id %in% node_ids, "darkorange", "steelblue")
  nodes_viz[, color.border := border_color_col] 
  
  
  # shape
  shape_col <- ifelse(nodes_viz$id %in% node_ids, "star", "dot")
  nodes_viz[, shape := shape_col] 
  
  # size of node according to centrality (here: degree)
  nodes_viz[, value := degree] 
  
  
  ### Edges
  
  # color
  color_col_edges <- ifelse(edges_viz$SAME_TAG == 1 & edges_viz$BIOGRID == 1, "blue", 
                            ifelse(edges_viz$SAME_TAG == 1, "darkgreen", "darkorange"))
  edges_viz$color <- color_col_edges
  
  # edges colour for candidates
  # color_col_edges <- ifelse(edges_viz$from %in% node_ids | edges_viz$to %in% node_ids, "orange", "lightskyblue")
  # edges_viz$color <- color_col_edges
  
  # label 
  label_col_edges <- ifelse(edges_viz$SAME_TAG == 1 & edges_viz$BIOGRID == 1, "SAME TAG + BIOGRID", 
                            ifelse(edges_viz$SAME_TAG == 1, "SAME TAG", "BIOGRID"))
  
  edges_viz$title <- label_col_edges
  
 
  
  ### Network output
  res_net <- visNetwork(nodes_viz, edges_viz, height = "500px") %>%
              visIgraphLayout() %>%
              visOptions(highlightNearest = list(enabled = TRUE, 
                                                 hover = TRUE), 
                         nodesIdSelection = TRUE)  
  
  
  res_list <- list("net" = res_net, 
                   "nodes" = nodes_viz, 
                   "edges" = edges_viz)
  
  return(res_list)
 
}


#' Function to show the gene-gene interaction network from rs_ids 
#'
#' @param rs_ids Vector with rs ids on which the network will be based on
#' @param dbSNP_dt Data table with the snp-gene mapping 
#' @param biogrid_dt Data table with gene-gene interaction information
#'
#' @return gene-gene interaction network based on biogrid 
#'
#' @examples
#' rs_ids <- c("rs41440544","rs41501252","rs16976638","rs16976644","rs10518828", "rs16976648", "rs767483060", "rs1989372")
#' show_gene_network_from_ids(rs_ids = rs_ids, dbSNP_dt = dbSNP_dt, biogrid_dt = biogrid_dt)
show_gene_network_from_ids <- function(rs_ids, dbSNP_dt, biogrid_dt){
  
  # Substitute rs in the rs_ids by ""
  rs_ids <- gsub("rs", "", rs_ids)
  
  # 1. SNPs to genes
  # rs ids to genes
  gene_nodes_dt <- create_rs_id_gene_table(rs_ids, dbSNP_dt)
  gene_nodes <- unique(unlist(strsplit(gene_nodes_dt$`gene/s`, ";")))
  
  
  # 2. Look in biogrid for edges 
  # get biogrid edges
  biogrid_edges_dt <- biogrid_dt[, .(`Official Symbol Interactor A`, `Official Symbol Interactor B`)]
  colnames(biogrid_edges_dt) <- c("from_gene", "to_gene")
  
  # get edges for genes
  gene_edges <- unique(biogrid_edges_dt[from_gene %in% gene_nodes & to_gene %in% gene_nodes])
  
  # gene node data table with ids 
  gene_nodes <- data.table(gene = str_sort(gene_nodes))
  gene_nodes[, id := seq(0, dim(gene_nodes)[1]-1)]
  
  
  # ids for edges
  from_gene_id <- sapply(gene_edges$from_gene, function(x) gene_nodes[gene == x, id])
  to_gene_id <- sapply(gene_edges$to_gene, function(x) gene_nodes[gene == x, id])
  gene_edges[, from := from_gene_id]
  gene_edges[, to := to_gene_id]
  
  
  ### Nodes
  colnames(gene_nodes)[1] <- "label"
  
  # label color
  font_color_col <- rep("black", dim(gene_nodes)[1])
  gene_nodes[, font.color := font_color_col]
  
  
  # shape
  shape_col <- rep("ellipse", dim(gene_nodes)[1])
  gene_nodes[, shape := shape_col]
  
  
  
  # with defaut layout
  net <- visNetwork(gene_nodes, gene_edges, height = "500px") %>%
          visIgraphLayout() %>%
          visNodes(scaling = list(min = 20, max = 50)) %>%
          visOptions(highlightNearest = list(enabled = T, hover = T),
                     nodesIdSelection = T)
  
  
  return(list("net" = net, 
              "gene_nodes_dt" = gene_nodes_dt))
}




#### Gene network function to keep
#####
#' Function to show the gene-gene interaction network from rs_ids 
#'
#' @param nodes Data table with columns: id, rs_id
#' @param edges Edges of the network
#' @param rs_ids Vector with rs ids on which the network will be based on
#' @param dbSNP_dt Data table with the snp-gene mapping 
#'
#' @return gene-gene interaction network 
#'
#' @examples 
#' rs_ids <- c("rs41440544","rs41501252","rs16976638","rs16976644","rs10518828",
#'             "rs16976648", "rs767483060", "rs1989372")
#' net_out <- show_gene_network_from_ids(nodes, edges, rs_ids, dbSNP_dt)

# show_gene_network_from_ids <- function(nodes, edges, rs_ids, dbSNP_dt){
#   
#     # Substitute rs in the rs_ids by ""
#     rs_ids <- gsub("rs", "", rs_ids)
#     
#     
#     # Get all ids for the given rs_ids
#     node_ids <- nodes[rs_id %in% rs_ids, id]
#     
#     
#     # edges from node ids
#     edges_viz <- edges[node1 %in% node_ids & node2 %in% node_ids][, .(node1, node2)]
#     colnames(edges_viz) <- c("from", "to")
#     
#     
#     # rs_ids from node_ids
#     from_rs <- c(sapply(edges_viz$from, function(x) nodes[id == x, rs_id]))
#     to_rs <- c(sapply(edges_viz$to, function(x) nodes[id == x, rs_id]))
#     
#     edges_viz[, from_rs := from_rs]
#     edges_viz[, to_rs := to_rs]
#     
#     
#     # gene/s from rs_ids
#     from_gene <- c(sapply(edges_viz$from_rs, function(x) paste(eval(get_gene_from_snp(rs_id = x, dbSNP_dt = dbSNP_dt)), collapse = ";")))
#     to_gene <- c(sapply(edges_viz$to_rs, function(x) paste(eval(get_gene_from_snp(rs_id = x, dbSNP_dt = dbSNP_dt)), collapse = ";")))
#     
#     edges_viz[, from_gene := from_gene]
#     edges_viz[, to_gene := to_gene]
#     
#     
#     # Gene edges only
#     edges_viz_gene <- edges_viz[, .(from_gene, to_gene)]
#     
#     
#     # Create all gene combinations when more than one gene are in an entry
#     gene_crosses <- data.table()
#     
#     if(length(from_gene) != 0){
#       
#       n <- length(from_gene)
#       
#       for (i in 1:n){
#         from_gene <- unlist(str_split(edges_viz_gene[i, c(from_gene)], ";"))
#         to_gene <- unlist(str_split(edges_viz_gene[i, c(to_gene)], ";"))
#         
#         tmp_cross <- crossing(from_gene, to_gene)
#         gene_crosses <- rbind(gene_crosses, tmp_cross)
#         
#       }
#     } 
#     
#     edges_viz_gene <- unique(gene_crosses)
#     
#     
#     ### Repeat process from rs_ids to genes with the rs_ids without edges
#     edges_viz_rs_ids <- unique(c(edges_viz$from_rs, edges_viz$to_rs))
#     wo_edges_rs_ids <- data.table(from_rs = rs_ids[!(rs_ids %in% edges_viz_rs_ids)], 
#                                   to_rs = rs_ids[!(rs_ids %in% edges_viz_rs_ids)])
#     
#     
#     # gene/s from rs_ids
#     from_gene <- c(sapply(wo_edges_rs_ids$from_rs, function(x) paste(eval(get_gene_from_snp(rs_id = x, dbSNP_dt = dbSNP_dt)), collapse = ";")))
#     to_gene <- c(sapply(wo_edges_rs_ids$to_rs, function(x) paste(eval(get_gene_from_snp(rs_id = x, dbSNP_dt = dbSNP_dt)), collapse = ";")))
#     
#     wo_edges_rs_ids[, from_gene := from_gene]
#     wo_edges_rs_ids[, to_gene := to_gene]
#     
#     # Gene edges only
#     wo_edges_gene <- wo_edges_rs_ids[, .(from_gene, to_gene)]
#     
#     # Create all gene combinations when more than one gene are in an entry
#     gene_crosses_wo <- data.table()
#     
#     if (length(from_gene) != 0){
#       
#       n <- length(from_gene)
#       
#       for (i in 1:n){
#         from_gene <- unlist(str_split(wo_edges_gene[i, c(from_gene)], ";"))
#         to_gene <- unlist(str_split(wo_edges_gene[i, c(to_gene)], ";"))
#         
#         tmp_cross <- crossing(from_gene, to_gene)
#         gene_crosses_wo <- rbind(gene_crosses_wo, tmp_cross)
#         
#       }
#     }
#     
#     wo_edges_gene <- unique(gene_crosses_wo)
#     
#     # check reflexive edges 
#     wo_edges_gene <- wo_edges_gene[! (from_gene == to_gene)]
#     
#     
#     # Combine all edges
#     edges_viz_gene <- rbind(edges_viz_gene, wo_edges_gene)
#     
#     
#     # nodes to visualize
#     nodes_viz_gene <- data.table(gene = str_sort(unique(c(edges_viz_gene$from_gene, edges_viz_gene$to_gene))))
#     nodes_viz_gene[, id := seq(0, dim(nodes_viz_gene)[1]-1)]
#     
#     # Ids for edges
#     from_gene_id <- sapply(edges_viz_gene$from_gene, function(x) nodes_viz_gene[gene == x, id])
#     to_gene_id <- sapply(edges_viz_gene$to_gene, function(x) nodes_viz_gene[gene == x, id])
#     edges_viz_gene[, from := from_gene_id]
#     edges_viz_gene[, to := to_gene_id]
#     
#     
#     
#     ### Nodes 
#     colnames(nodes_viz_gene)[1] <- "label"
#     
#     # label color 
#     font_color_col <- rep("black", dim(nodes_viz_gene)[1])
#     nodes_viz_gene[, font.color := font_color_col] 
#     
#     
#     # color
#     # color_col <- ifelse(nodes_viz_gene$id %in% node_ids, "orange", "lightskyblue")
#     # nodes_viz_gene[, color := color_col] 
#     
#     # shape
#     shape_col <- rep("ellipse", dim(nodes_viz_gene)[1])
#     nodes_viz_gene[, shape := shape_col]
#     
#     
#     
#     # with defaut layout
#     visNetwork(nodes_viz_gene, edges_viz_gene, height = "500px") %>%
#                   visIgraphLayout() %>%
#                   visNodes(scaling = list(min = 20, max = 50)) %>%
#                   visOptions(highlightNearest = list(enabled = T, hover = T), 
#                              nodesIdSelection = T) 
#     
#     
#     ## Table with rs_id and gene/s they belong to
#     # rs_id_gene_dt <- rbind(edges_viz[, -c("from", "to")], wo_edges_rs_ids)
#     # 
#     # rs_id_gene_dt <- unique(rbind(data.table(rs_id_gene_dt$from_rs, rs_id_gene_dt$from_gene), 
#     #                               data.table(rs_id_gene_dt$to_rs, rs_id_gene_dt$to_gene)))
#     # 
#     # colnames(rs_id_gene_dt) <- c("rs_id", "gene")
#     # 
#     # rs_id_gene_dt <- rs_id_gene_dt[order(gene)]
#     # 
#     
# }


