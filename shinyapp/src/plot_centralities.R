###############################################
#### Centrality plotting for given network ####
###############################################


# Function to turn adjacency list into rs_id degree data table 
adj_to_rs_degree <- function(adj_list){
  rs_id_degree_dt <- data.table(rs_id = as.integer(names(adj_list)), degree = sapply(adj_list, function(x) length(x)))
  return(rs_id_degree_dt)
}



#' Function to plot centrality distribution of whole SNP-SNP interaction network
#'
#' @param centr_dt Data table with rs_id 
#' @param centrality 
#' @param num_bins 
#' @param disease_name 
#' @param x_axis 
#'
#' @return Plot with distribution of given centrality 
#' @export
#'
#' @examples plot_centrality(nodes = nodes, centrality = "degree", num_bins = 50)
plot_centrality <- function(nodes, centrality = "", num_bins = 50){
  
  centr_plt <- ggplot() + 
    geom_histogram(data = nodes[! is.na(get(centrality)), .(rs_id, get(centrality))], 
                   aes(x = V2), bins = num_bins, fill = "grey")  +
    scale_y_log10() +
    theme_bw(base_size = 14) + 
    labs(x = eval(centrality), y = expression(log[10](count)))
  
  return(centr_plt)
}


#' Get first n candidate sets from result data table
#'
#' @param res_dt Data table with candidate information
#' @param n Number of entries to get back
#' @param wo_rs Boolean if rs should be remained or removed 
#'
#' @return Data table with list of rs_ids as elements 
#'
#' @examples
#' get_cands(res_dt = res_ext_dt)
get_cands <- function(res_dt, n = 10, wo_rs = TRUE){
  first_n_dt <- res_dt[1:n, .(RS_IDS)]
  if (wo_rs == TRUE){
    first_n_dt <- data.table(RS_IDS = gsub("rs", "", first_n_dt[, RS_IDS]))
  }
  first_n_dt <- sapply(first_n_dt, function(x) c(strsplit(x, split=";")))
  
  return(first_n_dt)
}




#' Create plot of whole distribution of a centrality across the 
#' SNP-SNP interaction network with marked input SNPs (usually candidate SNPs)
#'
#' @param nodes Data table with nodes and centrality information
#' @param rs_ids Rs_ids that should be marked
#' @param centrality Centrality to use for distribution
#' @param num_bins Number of bins
#'
#' @return Plot with marked SNPs in whole SNP-SNP interaction network 
#' distribution for given centrality 
#'
#' @examples
plot_centr_with_cand <- function(nodes, rs_ids, centrality, num_bins = 50){
  
  rs_id_centr_dt <- nodes[! is.na(get(centrality)), .(rs_id, get(centrality))]
  
  # add colnames for dt
  colnames(rs_id_centr_dt) <- c("rs_id", centrality)
  
  # remove rs from rs_ids 
  rs_ids <- gsub("rs", "", rs_ids)
  
  # candidate column
  cand_rs_id_centr_dt <- rs_id_centr_dt[rs_id %in% rs_ids]
  
  # mark candidates
  cand_rs_id_centr_dt$y <- rep(0, dim(cand_rs_id_centr_dt)[1])
  
  # Rename columns
  cand_rs_id_centr_dt$rs_id <- paste0("rs", cand_rs_id_centr_dt$rs_id)

  
  # create plot with marked candidates 
  cand_plt <- ggplot() + 
    geom_histogram(data = rs_id_centr_dt, aes(x = get(centrality)), bins = num_bins, fill = "grey") + 
    scale_y_log10() +
    geom_dotplot(data = cand_rs_id_centr_dt, aes(x=get(centrality), y=y, fill=rs_id), dotsize = (2/3), 
                 bins = num_bins, stackdir = "up", stackratio = 1.4, stackgroups=TRUE) +
    labs(x = eval(centrality), y = expression(log[10](count))) + 
    theme_bw(base_size = 14)
  
  return(cand_plt)
}




#' Function to get the degrees of the first n candidates and belonging to 
#' a candidate set 
#'
#' @param rs_id_degree_dt Data table with rs_id and degree
#' @param res_dt Data table with candidate information
#' @param n Number of entries to get back
#'
#' @return Data table with rs_id, degree and belonging to a candidate set 
#'
#' @examples get_cand_degrees(rs_id_degree_dt, res_dt, n = 10)
get_cand_degrees <- function(rs_id_degree_dt, res_dt, n = 10){
  cands <- get_cands(res_dt, n)
  rs_id_degree_dt <- rs_id_degree_dt[, .(rs_id, degree)]
  
  log_form <- ""
  for(i in 1:n) { 
    # candidate column
    tmp_cand_name <- paste0("cand_", i)
    tmp_cand <- as.numeric(cands[[i]])
    
    # filter for candidates in rs_ids of rs_id_degree_dt
    tmp_col <- rs_id_degree_dt$rs_id %in% tmp_cand
    rs_id_degree_dt[, paste0("cand_", i) := tmp_col]
    if(i != n){
      log_form <- paste0(log_form, tmp_cand_name, ",")
    } else {
      log_form <- paste0(log_form, tmp_cand_name, "")
    }
  }
  
  # create columns for n cands 
  cols <- unlist(strsplit(log_form, split=","))
  
  # filter for all candidate degrees and information in which candidate set they are 
  rs_id_degree_dt <- rs_id_degree_dt[rs_id_degree_dt[, do.call(pmax, .SD) > 0, .SDcols = cols]]
  
  return(rs_id_degree_dt[order(degree)])
}



# Function to plot centrality distribution of whole SSI with marked candidate sets
plot_centr_cands <- function(rs_id_degree_dt, res_dt, centrality = "", n = 10, num_bins = 300, 
                             disease_name = "", network_type = "", x_axis = "centrality", range_x = 10000){
  cands <- get_cands(res_dt, n)
  centr_cand_plts <- list()
  colnames(res_dt)[3] <- "SCORE"
  
  
  for(i in 1:n) { 
    # Candidate column
    tmp_score <- round(res_dt[i, SCORE], 2)
    tmp_cand_name <- paste0("cand_", i)
    tmp_cand <- as.numeric(cands[[i]])
    tmp_col <- rs_id_degree_dt$rs_id %in% tmp_cand
    rs_id_degree_dt[, paste0("cand_", i) := tmp_col]
    
    # Mark candidates
    mark_cand_x <- rs_id_degree_dt[get(tmp_cand_name) == T, get(centrality)]
    mark_cand_y <- rep(0, length((mark_cand_x)))
    mark_cand_dt <- data.table(x = mark_cand_x, y = mark_cand_y)
    
    
    # Degree candidate plots 
    # create plot
    tmp_plt <- ggplot() + 
      geom_histogram(data = rs_id_degree_dt, aes(x = get(centrality)), binwidth = num_bins, fill = "grey") + 
      ggtitle(paste("cand", i)) + 
      scale_y_log10() +
      geom_dotplot(data = mark_cand_dt, aes(x=x, y=y), dotsize = 2, fill = "blue", binwidth = num_bins) +
      labs(x = eval(x_axis), y = expression(log[10](count))) + 
      annotate("text", x = range_x, y = 1e5, label = paste("MLM:", tmp_score), hjust = 1, vjust = 1) + 
      theme_bw(base_size = 14)
    
    
    centr_cand_plts[[i]] <- tmp_plt
  }
  
  
  wrap_plots(centr_cand_plts) +
    plot_annotation(
      title = "DEGREE DISTRIBUTION WITH MARKED CANDIDATE SNP SETS",
      subtitle = paste0(network_type, ", " ,disease_name)
      # caption = ""
    )
  
}


