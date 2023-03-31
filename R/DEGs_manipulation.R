#' Perform Differential Expression Analysis
#' This function is used to perform DEA for single-cell RNA-seq data by QLF model in edgeR package including the cellular detection rate.
#' @param data A matrix of gene expression profile.
#' @param group Group(or cluster) assignment of every cells in columns of matrix.
#' @param method Which DE method will be used? Choices: edgeRQLF, edgeRQLFDetRate, MASTcpmDetRate, MASTtpmDetRate, MASTcpm, MASTtpm, limmatrend, limmavoom, ttest and wilcox.
#' @param species We only support `human` or `mouse` to choose the gene length data
#' from `simutils::hs_gene_length` or `simutils::mm_gene_length` and normalize the data
#' by TPM. Users must set this parameter when using `MASTtpmDetRate`, `MASTtpm`, `ttest` or `wilcox` method.
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest topTags cpm
#' @importFrom limma topTable lmFit eBayes voom
#' @importFrom stats model.matrix p.adjust t.test
#' @importFrom methods new
#' @return A dataframe of DEA result.
#' @export
#' @references
#' Soneson C, Robinson M D. Bias, robustness and scalability in single-cell differential expression analysis. Nature methods, 2018, 15(4): 255-261. <https://doi.org/10.1038/nmeth.4612>
#' Github URL: <https://github.com/csoneson/conquer_comparison>
perform_DEA <- function(
  data,
  group,
  method,
  species = NULL
){
  ## Filter
  filter_index <- rowSums(data) > 0
  message(paste0(nrow(data) - sum(filter_index), " genes are removed when filtering"))
  data <- data[rowSums(data) > 0, ]
  ## paired groups
  group_unique <- unique(group)
  group_paired <- utils::combn(group_unique, 2)
  ## blank result list
  result_list <- list()
  for(i in 1:ncol(group_paired)){
    group_candidate <- group_paired[, i]
    index <- group %in% group_candidate
    sub_group <- group[index]
    sub_data <- data[, index]
    message("--------------------------------------------------")
    message(glue::glue("Performing DEA with the {i}/{ncol(group_paired)} paired of groups."))

    if(method == "edgeRQLF"){
      message("Performing DEA using edgeR QLF model")
      timing <- system.time({
        dge <- edgeR::DGEList(sub_data, group = sub_group)
        dge <- edgeR::calcNormFactors(dge)
        design <- stats::model.matrix(~ sub_group)
        dge <- edgeR::estimateDisp(dge, design = design)
        fit <- edgeR::glmQLFit(dge, design = design)
        qlf <- edgeR::glmQLFTest(fit)
        tt <- edgeR::topTags(qlf, n = Inf)
      })
      result_list[[paste0(group_candidate, collapse = "vs")]] <- tt$table
    }

    if(method == "edgeRQLFDetRate"){
      message("Performing DEA using edgeR QLF model including the cellular detection rate")
      timing <- system.time({
        dge <- edgeR::DGEList(sub_data, group = sub_group)
        dge <- edgeR::calcNormFactors(dge)
        cdr <- scale(colMeans(sub_data > 0))
        design <- stats::model.matrix(~ cdr + sub_group)
        dge <- edgeR::estimateDisp(dge, design = design)
        fit <- edgeR::glmQLFit(dge, design = design)
        qlf <- edgeR::glmQLFTest(fit)
        tt <- edgeR::topTags(qlf, n = Inf)
      })
      result_list[[paste0(group_candidate, collapse = "vs")]] <- tt$table
    }

    if(method == "MASTcpmDetRate"){
      if(!requireNamespace("MAST", quietly = TRUE)){
        stop("Package \"MAST\" must be installed by \"BiocManager::install('MAST')\".")
      }
      message("Performing DEA using MAST including the cellular detection rate with CPM data")
      timing <- system.time({
        names(sub_group) <- colnames(sub_data)
        cdr <- scale(colMeans(sub_data > 0))
        dge <- edgeR::DGEList(counts = sub_data)
        dge <- edgeR::calcNormFactors(dge)
        cpms <- edgeR::cpm(dge)
        sca <- MAST::FromMatrix(exprsArray = log2(cpms + 1),
                                cData = data.frame(wellKey = names(sub_group),
                                                   group = sub_group, cdr = cdr))
        zlmdata <- MAST::zlm(~ cdr + group, sca)
        mast <- MAST::lrTest(zlmdata, "group")
        df <- data.frame(PValue = mast[, "hurdle", "Pr(>Chisq)"],
                         FDR = stats::p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH"),
                         row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
      })
      result_list[[paste0(group_candidate, collapse = "vs")]] <- df
    }

    if(method == "MASTtpmDetRate"){
      if(!requireNamespace("MAST", quietly = TRUE)){
        stop("Package \"MAST\" must be installed by \"BiocManager::install('MAST')\".")
      }
      message("Performing DEA using MAST including the cellular detection rate with TPM data")
      tpm <- simutils::normalization_simutils(data = sub_data,
                                              norm_method = "TPM",
                                              species = species)
      timing <- system.time({
        names(sub_group) <- colnames(tpm)
        cdr <- scale(colMeans(tpm > 0))
        sca <- MAST::FromMatrix(exprsArray = log2(tpm + 1),
                                cData = data.frame(wellKey = names(sub_group),
                                                   group = sub_group, cdr = cdr))
        zlmdata <- MAST::zlm(~ cdr + group, sca)
        mast <- MAST::lrTest(zlmdata, "group")
        df <- data.frame(PValue = mast[, "hurdle", "Pr(>Chisq)"],
                         FDR = stats::p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH"),
                         row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
      })
      result_list[[paste0(group_candidate, collapse = "vs")]] <- df
    }

    if(method == "MASTcpm"){
      if(!requireNamespace("MAST", quietly = TRUE)){
        stop("Package \"MAST\" must be installed by \"BiocManager::install('MAST')\".")
      }
      message("Performing DEA using MAST with CPM data")
      timing <- system.time({
        names(sub_group) <- colnames(sub_data)
        dge <- edgeR::DGEList(counts = sub_data)
        dge <- edgeR::calcNormFactors(dge)
        cpms <- edgeR::cpm(dge)
        sca <- MAST::FromMatrix(exprsArray = log2(cpms + 1),
                                cData = data.frame(wellKey = names(sub_group),
                                                   group = sub_group))
        zlmdata <- MAST::zlm(~ group, sca)
        mast <- MAST::lrTest(zlmdata, "group")
        df <- data.frame(PValue = mast[, "hurdle", "Pr(>Chisq)"],
                         FDR = stats::p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH"),
                         row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
      })
      result_list[[paste0(group_candidate, collapse = "vs")]] <- df
    }

    if(method == "MASTtpm"){
      if(!requireNamespace("MAST", quietly = TRUE)){
        stop("Package \"MAST\" must be installed by \"BiocManager::install('MAST')\".")
      }
      message("Performing DEA using MAST with TPM data")
      tpm <- simutils::normalization_simutils(data = sub_data,
                                              norm_method = "TPM",
                                              species = species)
      timing <- system.time({
        names(sub_group) <- colnames(tpm)
        sca <- MAST::FromMatrix(exprsArray = log2(tpm + 1),
                                cData = data.frame(wellKey = names(sub_group),
                                                   group = sub_group))
        zlmdata <- MAST::zlm(~ group, sca)
        mast <- MAST::lrTest(zlmdata, "group")
        df <- data.frame(PValue = mast[, "hurdle", "Pr(>Chisq)"],
                         FDR = stats::p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH"),
                         row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
      })
      result_list[[paste0(group_candidate, collapse = "vs")]] <- df
    }

    if(method == "limmatrend"){
      message("Performing DEA using limma-trend")
      timing <- system.time({
        dge <- edgeR::DGEList(sub_data, group = sub_group)
        dge <- edgeR::calcNormFactors(dge)
        design <- stats::model.matrix(~ sub_group)
        y <- methods::new("EList")
        y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
        fit <- limma::lmFit(y, design = design)
        fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
        tt <- limma::topTable(fit, n = Inf, adjust.method = "BH")
        colnames(tt)[c(4, 5)] <- c("PValue", "FDR")
      })
      result_list[[paste0(group_candidate, collapse = "vs")]] <- tt
    }

    if(method == "limmavoom"){
      message("Performing DEA using limma-voom")
      timing <- system.time({
        dge <- edgeR::DGEList(sub_data, group = sub_group)
        dge <- edgeR::calcNormFactors(dge)
        design <- model.matrix(~ sub_group)
        vm <- limma::voom(dge, design = design, plot = FALSE)
        fit <- limma::lmFit(vm, design = design)
        fit <- limma::eBayes(fit)
        tt <- limma::topTable(fit, n = Inf, adjust.method = "BH")
        colnames(tt)[c(4, 5)] <- c("PValue", "FDR")
      })
      result_list[[paste0(group_candidate, collapse = "vs")]] <- tt
    }

    if(method == "ttest"){
      message("Performing DEA using t-test")
      tpm <- simutils::normalization_simutils(data = sub_data,
                                              norm_method = "TPM",
                                              species = species)
      timing <- system.time({
        tmm <- edgeR::calcNormFactors(tpm)
        tpmtmm <- edgeR::cpm(tpm, lib.size = tmm * colSums(tpm))
        logtpm <- log2(tpmtmm + 1)
        idx <- seq_len(nrow(logtpm))
        names(idx) <- rownames(logtpm)
        ttest_p <- sapply(idx, function(i) {
          stats::t.test(logtpm[i, ] ~ sub_group)$p.value
        })
      })
      df <- data.frame(PValue = ttest_p,
                       FDR = stats::p.adjust(ttest_p, method = "BH"),
                       row.names = names(ttest_p))
      result_list[[paste0(group_candidate, collapse = "vs")]] <- df
    }

    if(method == "wilcox"){
      message("Performing DEA using wilcox-test")
      tpm <- simutils::normalization_simutils(data = sub_data,
                                              norm_method = "TPM",
                                              species = species)
      timing <- system.time({
        tmm <- edgeR::calcNormFactors(tpm)
        tpmtmm <- edgeR::cpm(tpm, lib.size = tmm * colSums(tpm))
        idx <- 1:nrow(tpmtmm)
        names(idx) <- rownames(tpmtmm)
        wilcox_p <- sapply(idx, function(i) {
          stats::wilcox.test(tpmtmm[i, ] ~ sub_group)$p.value
        })
      })
      df <- data.frame(PValue = wilcox_p,
                       FDR = stats::p.adjust(wilcox_p, method = "BH"),
                       row.names = names(wilcox_p))
      result_list[[paste0(group_candidate, collapse = "vs")]] <- df
    }
  }
  return(result_list)
}



#' Proportion of Matched DEGs
#'
#' @param sim_data The simulated dataset to be evaluated.
#' @param group Group(or cluster) assignment of every cells in columns of matrix.
#' @param group_combn The combinations of two groups generated by [utils::combn()] function.
#' @param sim_DEGs A list with `xxxvsxxx` format of names containing the DEGs of different paired groups.
#' @param DEA_method The method used for differential expression analysis. edgeRQLF, edgeRQLFDetRate, MASTcpmDetRate,
#' MASTcpm, limmatrend and limmavoom.
#' @importFrom purrr map
#' @importFrom BiocGenerics Reduce intersect
#' @export
true_DEGs_proportion <- function(
    sim_data,
    group,
    group_combn,
    sim_DEGs,
    DEA_method
){
  true_prop <- list()
  DEGs_num <- c()
  DEA <- list()
  DEGs_total <- length(unique(unlist(sim_DEGs)))
  for(i in 1:ncol(group_combn)){
    group_compare <- group_combn[, i]
    compare_name <- paste0(group_compare, collapse = "vs")
    sub_data <- sim_data[, which(group %in% group_compare)]
    sub_group <- group[which(group %in% group_compare)]
    sub_DEGs <- sim_DEGs[[compare_name]]
    all_method_DEGs <- purrr::map(DEA_method, function(method){
      DEA_result <- simutils::perform_DEA(data = sub_data,
                                          group = sub_group,
                                          method = method)
      list("de_genes" = rownames(DEA_result[[1]])[DEA_result[[1]]$"PValue" < 0.05],
           "DEA_result" = DEA_result[[1]])
    }) %>% setNames(DEA_method)
    if(length(DEA_method) > 1){
      intersect_genes_methods <- list()
      for(method in DEA_method){
        intersect_genes_methods[[method]] <- all_method_DEGs[[method]][["de_genes"]]
      }
      intersect_genes_methods <- BiocGenerics::Reduce(x = all_method_DEGs, f = intersect)
    }else{
      intersect_genes_methods <- all_method_DEGs[[1]][["de_genes"]]
    }
    intertect_genes <- BiocGenerics::intersect(intersect_genes_methods, sub_DEGs)
    prop <- length(intertect_genes)/length(sub_DEGs)
    DEGs_num <- append(DEGs_num, length(sub_DEGs))
    true_prop[compare_name] <- prop
    DEA[[compare_name]] <- all_method_DEGs
  }
  weighted_true_prop <- sum(unname(unlist(true_prop)) * (DEGs_num / sum(DEGs_num)))
  return_list <- dplyr::lst(true_prop,
                            DEA,
                            DEGs_num,
                            DEGs_total,
                            weighted_true_prop)
  return(return_list)
}



#' Test if the p-value is Uniformly Distributed
#'
#' @param PValue The p-value vector.
#' @export
#' @examples
#' pvalue <- runif(1000)
#' result <- test_uni_distribution(pvalue)
test_uni_distribution <- function(
    PValue
){
  if(!requireNamespace("spgs", quietly = TRUE)){
    message("spgs is not installed on your device...")
    message("Installing spgs...")
    install.packages("spgs")
  }
  result <- spgs::chisq.unif.test(PValue, bins = 20)
  if(result$p.value < 0.05){
    score <- 0
  }else{
    score <- 1
  }
  return(dplyr::lst(result,
                    score))
}



#' Summarize the Ability of Simulating DEGs
#'
#' @param count_matrix. A matrix or a list with names of gene expression data.
#' @param group. A vector of characters which indicate which group a cell belongs to. A list is also supported when many matrices of data is input in `count_matrix` parameter.
#' @param DEGs. A list of DEGs with the names of `xxxvsxxx`. Note that the names of DEGs are in the rownames of the matrix or the dataframe and the names of `xxx` is in the `group` characters. If you have input the lists of `group` and `count_matrix`, the every sub-list in DEGs list should match with every data in `count_matrix` and `group` at the same index or position.
#' @importFrom utils combn
#' @importFrom stringr str_split
#' @importFrom dplyr lst
#' @importFrom stats setNames
#' @importFrom purrr map
#' @return A list
#' @export
#'
calculate_DEGs_properties <- function(
  count_matrix,
  group,
  DEGs
){
  ### Check input
  if(is.matrix(count_matrix)){
    count_matrix <- list(data = count_matrix)
  }

  if(!is.list(group)){
    group <- list(group = group)
  }
  if(length(count_matrix) != length(group)){
    stop("The length of input data is not equal to the length of correponding group information")
  }
  for(id in 1:length(count_matrix)){
    if(ncol(count_matrix[[id]]) != length(group[[id]])){
      stop(paste0("The cell number is not equal to the length of group information in the ", id, "th dataset"))
    }
  }

  if(length(count_matrix) == 1 & length(DEGs) == 3){
    if(all(stringr::str_detect(names(DEGs), pattern = "vs"))){
      DEGs <- list(DEGs = DEGs)
    }else{
      stop("The names of list DEGs are not valid or the length of data you have input is nout equal to that of DEGs list")
    }
  }
  data_length <- length(count_matrix)
  data_names <- names(count_matrix)
  DEGs_evaluation <- purrr::map(
    .x = 1:data_length,
    .f = function(data_id){
      data <- count_matrix[[data_id]]
      if(is.data.frame(data)){
        data <- as.matrix(data)
      }
      group <- group[[data_id]]
      DEGs <- DEGs[[data_id]]
      group_combn <- utils::combn(unique(group), 2)

      ### Distribution of P values
      distribution_score <- c()
      valid_DEGs_distribution <- list()
      for(i in 1:length(DEGs)){
        sub_DGEs <- DEGs[[i]]
        compare_name <- names(DEGs)[i]
        group1 <- stringr::str_split(compare_name,
                                     pattern = "vs",
                                     simplify = TRUE)[1]
        group2 <- stringr::str_split(compare_name,
                                     pattern = "vs",
                                     simplify = TRUE)[2]

        ### Distribution
        message("Distribution of null data...")
        index_DEGs <- which(rownames(data) %in% sub_DGEs)
        index1 <- which(group %in% group1)
        index2 <- which(group %in% group2)
        sub_data <- data[-index_DEGs, c(index1, index2)]
        sub_group <- c(rep(group1, length(index1)), rep(group2, length(index2)))
        sub_DEA_result <- perform_DEA(data = sub_data,
                                      group = sub_group,
                                      method = "edgeRQLFDetRate")
        p_values <- sub_DEA_result[[1]][["PValue"]]
        uniform_result <- test_uni_distribution(p_values)
        distribution_score <- append(distribution_score, uniform_result[["score"]])
        valid_DEGs_distribution[[compare_name]] <- sub_DEA_result[[1]]
      }

      ### True proportions of DEGs
      message("True proportions of DEGs...")
      DEGs_result <- true_DEGs_proportion(sim_data = data,
                                          group = group,
                                          group_combn = utils::combn(unique(group), 2),
                                          sim_DEGs = DEGs,
                                          DEA_method = "edgeRQLFDetRate")
      true_proportion <- DEGs_result[["weighted_true_prop"]]

      ### SVM
      de_genes <- unique(unlist(DEGs))
      SVM_result <- model_predict(data = data,
                                  group = group,
                                  de_genes = DEGs,
                                  method = "SVM")
      AUC <- as.numeric(SVM_result$roc$auc)
      Accuracy <- unname(SVM_result[["conf_matrix"]][["overall"]][1])
      if(length(unique(group)) == 2){
        Precision <- unname(SVM_result[["conf_matrix"]][["byClass"]]["Precision"])
        Recall <- unname(SVM_result[["conf_matrix"]][["byClass"]]["Recall"])
        F1 <- unname(SVM_result[["conf_matrix"]][["byClass"]]["F1"])
      }else{
        Precision <- mean(SVM_result[["conf_matrix"]][["byClass"]][, "Precision"], na.rm = TRUE)
        Recall <- mean(SVM_result[["conf_matrix"]][["byClass"]][, "Recall"], na.rm = TRUE)
        F1 <- mean(SVM_result[["conf_matrix"]][["byClass"]][, "F1"], na.rm = TRUE)
      }

      dplyr::lst(
        distribution_score = list(DEA_result_combinations = valid_DEGs_distribution,
                                  distribution_score = distribution_score),
        true_proportions = list(DEA_result_combinations = DEGs_result,
                                true_proportion = true_proportion),
        SVM_result <- list(SVM_result = SVM_result,
                           AUC = AUC,
                           Accuracy = Accuracy,
                           Precision = Precision,
                           Recall = Recall,
                           F1 = F1)
      )
    }
  ) %>% stats::setNames(data_names)
}





