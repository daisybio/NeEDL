#####################
#### Plot Scores ####
#####################


#' Plot the search score over time 
#'
#' @param search_score_dt Data table with search score information
#'
#' @return Plot with search scores starting from the second value as the 
#' first value is too high to plot
#' 
#' @examples plot_search_score(search_score_dt)
plot_search_score <- function(search_score_dt){
  
  colnames(search_score_dt) <- c("time_in_ms", "score")
  
  search_score_plt <- ggplot(search_score_dt, aes(x = time_in_ms, y = score)) + 
                          geom_step() + 
                          labs(x = "time (ms)", y = "MLM") +
                          theme_bw(base_size=14)
  
  return(search_score_plt)
}


#' Plot all scores e.g. bayesian, penetrance, regression, ...
#'
#' @param res_dt  Data table with candidate SNP sets with rs_ids and scores
#' @param score_names score names like bayesian, penetrance, regression for categories
#'
#' @return Boxplot for each model
#' @export
#'
#' @examples plot_set_scores(res_dt, score_names)
plot_set_scores <- function(res_dt, score_names){
  
  # Get all columns with scores defined by score_names
  score_dt <- res_dt[, colnames(res_dt) %in% score_names, with=FALSE]
  
  # Melt score_dt for plotting with facet_wrap
  score_dt <- melt(data = score_dt, 
                   measure.vars = colnames(score_dt), 
                   variable.name = "model_score", value.name = "score")
  
  # Plot boxplots using ggplot
  scores_plt <- ggplot(data = score_dt, aes(x = factor(0), y = score)) + 
                    geom_boxplot() + 
                    theme_bw(base_size = 14) + 
                    theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank()) + 
                    facet_wrap(~model_score, scales = "free")
  
  return(scores_plt)
}


#' Plot statistical scores of the individual candidate SNPs
#'
#' @param rs_ids rs_ids for which the scores should be shown in the plot
#' @param snp_scores_dt data table with score information
#' @param score_name name of the score as string 
#' @param show_gene boolean to display gene 
#'
#' @return Plot with rs_ids on the x axis and corresponding score on the y axis
#'
#' @examples
#' names(snp_scores_dt)[names(snp_scores_dt) == 'RS_IDS'] <- "rs_id"
#' rs_ids <- c("rs7596121", "rs7589521", "rs16841838", "rs41501252", "rs4947984")
#' score_name <- "REGRESSION_NLL"
#' plot_snp_scores(rs_ids, snp_scores_dt, score_name)
plot_snp_scores <- function(rs_ids, snp_scores_dt, score_name, show_gene = FALSE){
  # remove rs from rs ids 
  rs_ids <- gsub("rs", "", rs_ids)
  
  # get all entries for given rs ids 
  snp_scores_dt <- snp_scores_dt[rs_id %in% rs_ids]
 
  # add rs id 
  snp_scores_dt$rs_id <- paste0("rs", snp_scores_dt$rs_id)
 
  # Rename columns
  names(snp_scores_dt)[names(snp_scores_dt) == eval(score_name)] <- "score"
  
  # plot scores for rsid
  if(show_gene == FALSE){
    scores_plt <- ggplot(data = snp_scores_dt[order(score)], 
                         aes(x = rs_id, 
                             y = score, 
                             text = paste0(
                               "rs_id: ", rs_id, "<br>",
                               eval(score_name), ": ", score, "<br>"
                             ))) + 
      geom_point(size = 2, color = "#2596be") + 
      labs(x = "rs_id", y = score_name) +
      scale_x_discrete(limits=snp_scores_dt$rs_id) + 
      # scale_y_continuous(expand = c(0, 0), limits = c(0, max(snp_scores_dt[, score]) + 500)) + 
      theme_bw(base_size = 14) +  
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    
    return(ggplotly(scores_plt, tooltip = c("text")))
    
  } else {
    scores_plt <- ggplot(data = snp_scores_dt[order(score)], 
                         aes(x = rs_id, 
                             y = score, 
                             text = paste0(
                               "rs_id: ", rs_id, "<br>",
                               eval(score_name), ": ", score, "<br>", 
                               "gene/s: ", `gene/s`
                             ))) + 
      geom_point(size = 2, color = "#2596be") + 
      labs(x = "rs_id", y = score_name) +
      scale_x_discrete(limits=snp_scores_dt$rs_id) + 
      # scale_y_continuous(expand = c(0, 0), limits = c(0, max(snp_scores_dt[, score]) + 500)) + 
      theme_bw(base_size = 14) +  
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    
    return(ggplotly(scores_plt, tooltip = c("text")))
  }
 
  
}


#' Plot statistical scores of the individual candidate SNPs
#'
#' @param rs_ids rs_ids for which the scores should be shown in the plot
#' @param snp_scores_dt data table with score information
#' @param score_name name of the score as string 
#'
#' @return Plot with rs_ids on the x axis and corresponding score on the y axis
#'
#' @examples
#' names(snp_scores_dt)[names(snp_scores_dt) == 'RS_IDS'] <- "rs_id"
#' rs_ids <- c("rs7596121", "rs7589521", "rs16841838", "rs41501252", "rs4947984")
#' score_name <- "REGRESSION_NLL"
#' plot_snp_scores(rs_ids, snp_scores_dt, score_name)
plot_snp_scores_static <- function(rs_ids, snp_scores_dt, score_name){
  # remove rs from rs ids 
  rs_ids <- gsub("rs", "", rs_ids)
  
  # get all entries for given rs ids 
  snp_scores_dt <- snp_scores_dt[rs_id %in% rs_ids]
  
  # add rs id 
  snp_scores_dt$rs_id <- paste0("rs", snp_scores_dt$rs_id)
  
  # Rename columns
  names(snp_scores_dt)[names(snp_scores_dt) == eval(score_name)] <- "score"
  
  # plot scores for rsid
  scores_plt <- ggplot(data = snp_scores_dt[order(score)], 
                       aes(x = rs_id, 
                           y = score, 
                           text = paste0(
                             "rs_id: ", rs_id, "<br>",
                             eval(score_name), ": ", score, "<br>"
                           ))) + 
    geom_point(size = 3, color = "#2596be") + 
    labs(x = "rs_id", y = score_name) +
    scale_x_discrete(limits=snp_scores_dt$rs_id) + 
    # scale_y_continuous(expand = c(0, 0), limits = c(0, max(snp_scores_dt[, score]) + 500)) + 
    theme_bw(base_size = 14) +  
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  return(scores_plt)
   
  
  
}


