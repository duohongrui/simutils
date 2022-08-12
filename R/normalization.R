#' FPKM or TPM Normalization Method in simutils
#'
#' @param data A matrix with cells on columns and genes on rows.
#' @param norm_method FPKM or TPM
#' @param gene_length A numeric vector of gene lengths. If NULL, the gene length data will be imported by [simutils::hs_gene_length] or [simutils::mm_gene_length] function and users must set `species` parameters below.
#' @param species `human` or `mouse`.
#' @export
normalization_simutils <- function(
  data,
  norm_method = "FPKM",
  gene_length = NULL,
  species = "human"
){
  ## gene length data
  if(is.null(gene_length)){
    if(species == "human"){
      gene <- simutils::hs_gene_length
    }
    if(species == "mouse"){
      gene <- simutils::mm_gene_length
    }
  }
  ## gene form
  gene_name <- rownames(data)[1:10]
  if(all(startsWith(gene_name, "ENS"))){
    gene_formal <- "ensembl_id"
  }else{
    if(all(grepl(pattern = "^[1-9]", gene_name))){
      gene_formal <- "gene_id"
    }else{
      gene_formal <- "symbol"
    }
  }
  ## filter miss-match gene
  intersect_gene <- intersect(rownames(data), gene[, gene_formal])
  message(paste0(nrow(data)-length(intersect_gene), " genes are removed due to the miss-match."))
  data <- data[intersect_gene, ]
  ## FPKM
  index <- match(rownames(data), gene[, gene_formal])
  gene_length <- gene$length[index]
  data_by_length <- data/gene_length
  if(norm_method == "FPKM"){
    norm_data <- t(t(data_by_length) / colSums(data)) * 10^9
  }
  if(norm_method == "TPM"){
    norm_data <- t(t(data_by_length) / colSums(data_by_length)) * 10^6
  }
  return(norm_data)
}

