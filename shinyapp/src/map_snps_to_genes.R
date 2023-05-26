###########################
#### Map SNPs to genes ####
###########################


#' Function to get gene to which the given SNP corresponds 
#'
#' @param rs_id rs_id of the SNP for which the gene should be found
#' @param dbSNP_dt data table from dbSNP with SNP/gene mapping information
#'
#' @return gene for the given input SNP
#' @example get_gene_from_snp("548422872", dbSNP_dt)

get_gene_from_snp <- function(rs_id, dbSNP_dt){
  # Replace rs in rs_ids by ""
  rs_id <- as.numeric(gsub("rs", "", rs_id))
  
  # Create regex using the three required cases
  rs_id_mod_1 <- paste0("^", rs_id, ";")
  rs_id_mod_2 <- paste0(";", rs_id, ";")
  rs_id_mod_3 <- paste0(";", rs_id, "$")
  reg_expr <- paste0(rs_id_mod_1, "|" ,rs_id_mod_2, "|", rs_id_mod_3)
  
  # Look for gene/genes using the created regex
  gene <- dbSNP_dt[grepl(reg_expr, dbSNP_dt$SNPs), gene]
  
  return(gene)
}



#' Function to create a genes to SNP mapping from a SNPs to gene mapping
#'
#' @param dbSNP_dt: dbSNP with gene and corresponding SNP list per row
#'
#' @return genes to SNP mapping
#'
#' @examples
#' dbSNP_dt_snp_gene <- create_snp_to_gene_mapping(dbSNP_dt)
#' fwrite(dbSNP_dt_snp_gene, "/nfs/proj/GenEpiSeeker/genepiseeker_dev/data/dbSNP/inc_pseudogenes/snps_restruc_full_inc_pseudo_genes_to_snp.csv")

create_snp_to_gene_mapping <- function(dbSNP_dt){
  
  snp_to_gene_dt <- rbindlist(lapply(dbSNP_dt$gene, function(tmp_gene){
    
    # get row for given gene
    tmp_dbSNP_dt <- dbSNP_dt[gene == tmp_gene]
    
    # get snp list for given gene
    tmp_SNPs <- strsplit(as.character(tmp_dbSNP_dt$SNPs), split = ";", fixed = TRUE)
    
    # create data table with snp list and given gene (reused)
    tmp_SNPs_gene_dt <- data.table("SNP" = unlist(tmp_SNPs), 
                                   "genes" = tmp_gene)
    
    
  }))
  
  
  # group table using SNP as id and collapse corresponding genes 
  snp_to_gene_dt <- snp_to_gene_dt[, paste(genes, collapse=";"), by = SNP] 
  
  colnames(snp_to_gene_dt) <- c("SNP", "genes")
  
  return(snp_to_gene_dt)
  
}




#' Create rs_id gene table using dbSNP 
#'
#' @param rs_ids Vector with rs_ids from which the genes should be found
#' @param dbSNP_dt Lookup table with snp-gene mapping
#'
#' @return Data table with rs_id and gene/s
#'
#' @examples
#' rs_ids <- c("rs41440544","rs41501252","rs16976638","rs16976644","rs10518828","rs16976648")
#' create_rs_id_gene_table(rs_ids, dbSNP_dt)
#

create_rs_id_gene_table <- function(rs_ids, dbSNP_dt){
  # Substitute rs in the rs_ids by ""
  rs_ids <- gsub("rs", "", rs_ids)
  
  # Get genes to rs_ids
  genes_from_rs_ids <- c(sapply(rs_ids, function(x) paste(eval(get_gene_from_snp(rs_id = x, dbSNP_dt = dbSNP_dt)), collapse = ";")))
  
  # Create results table
  rs_id_gene_dt <- data.table("rs_id" = rs_ids, 
                              "gene/s" = genes_from_rs_ids)
  
  # Sort result table by gene lexicographically 
  rs_id_gene_dt <- rs_id_gene_dt[order(`gene/s`)]
  
  
  return(rs_id_gene_dt)
}






#' Create rs_id gene centrality table using dbSNP
#'
#' @param rs_id_centr_dt Data table with rs_id and centralities 
#' @param rs_ids Vector with rs_ids from which the genes should be found
#' @param dbSNP_dt Lookup table with snp-gene mapping
#'
#' @return Data table with rs_id, gene, centralities
#'
#' @examples
#' rs_ids <- c("rs41440544","rs41501252","rs16976638","rs16976644","rs10518828","rs16976648")
#' create_rs_id_gene_centr_table(nodes, rs_ids, dbSNP_dt)

create_rs_id_gene_centr_table <- function(rs_id_centr_dt, rs_ids, dbSNP_dt){

  # Substitute rs in the rs_ids by ""
  rs_ids <- gsub("rs", "", rs_ids)
  
  # Create rs_id gene table
  rs_id_gene_dt <- create_rs_id_gene_table(rs_ids, dbSNP_dt)
  
  # Create rs_id centrality table with only the rs_ids from above
  rs_id_gene_centr_dt <- rs_id_centr_dt[rs_id %in% rs_ids]
  
  # Remove id column
  rs_id_gene_centr_dt[, id := NULL]
  
  # Merge rs_id gene table with rs_id centrality table
  rs_id_gene_centr_dt <- merge(rs_id_gene_dt, rs_id_gene_centr_dt, by ="rs_id")
  
  rs_id_gene_centr_dt$rs_id <- sapply(rs_id_gene_centr_dt$rs_id, function(x) paste0("rs", x))
  
  return(rs_id_gene_centr_dt)
}

