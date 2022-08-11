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
#' @importFrom MAST FromMatrix zlm lrTest
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
  group <- as.numeric(as.factor(group))
  if(method == "edgeRQLF"){
    message("Performing DEA using edgeR QLF model")
    timing <- system.time({
      dge <- edgeR::DGEList(data, group = group)
      dge <- edgeR::calcNormFactors(dge)
      design <- stats::model.matrix(~ group)
      dge <- edgeR::estimateDisp(dge, design = design)
      fit <- edgeR::glmQLFit(dge, design = design)
      qlf <- edgeR::glmQLFTest(fit)
      tt <- edgeR::topTags(qlf, n = Inf)
    })
    return(tt$table)
  }

  if(method == "edgeRQLFDetRate"){
    message("Performing DEA using edgeR QLF model including the cellular detection rate")
    timing <- system.time({
      dge <- edgeR::DGEList(data, group = group)
      dge <- edgeR::calcNormFactors(dge)
      cdr <- scale(colMeans(data > 0))
      design <- stats::model.matrix(~ cdr + group)
      dge <- edgeR::estimateDisp(dge, design = design)
      fit <- edgeR::glmQLFit(dge, design = design)
      qlf <- edgeR::glmQLFTest(fit)
      tt <- edgeR::topTags(qlf, n = Inf)
    })
    return(tt$table)
  }

  if(method == "MASTcpmDetRate"){
    message("Performing DEA using MAST including the cellular detection rate with CPM data")
    timing <- system.time({
      names(group) <- colnames(data)
      cdr <- scale(colMeans(data > 0))
      dge <- edgeR::DGEList(counts = data)
      dge <- edgeR::calcNormFactors(dge)
      cpms <- edgeR::cpm(dge)
      sca <- MAST::FromMatrix(exprsArray = log2(cpms + 1),
                              cData = data.frame(wellKey = names(group),
                                                 group = group, cdr = cdr))
      zlmdata <- MAST::zlm(~ cdr + group, sca)
      mast <- MAST::lrTest(zlmdata, "group")
      df <- data.frame(PValue = mast[, "hurdle", "Pr(>Chisq)"],
                       FDR = stats::p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH"),
                       row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
    })
    return(df)
  }

  if(method == "MASTtpmDetRate"){
    message("Performing DEA using MAST including the cellular detection rate with TPM data")
    tpm <- simutils::normalization_simutils(data = data,
                                            norm_method = "TPM",
                                            species = species)
    timing <- system.time({
      names(group) <- colnames(tpm)
      cdr <- scale(colMeans(tpm > 0))
      sca <- MAST::FromMatrix(exprsArray = log2(tpm + 1),
                              cData = data.frame(wellKey = names(group),
                                                 group = group, cdr = cdr))
      zlmdata <- MAST::zlm(~ cdr + group, sca)
      mast <- MAST::lrTest(zlmdata, "group")
      df <- data.frame(PValue = mast[, "hurdle", "Pr(>Chisq)"],
                       FDR = stats::p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH"),
                       row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
    })
    return(df)
  }

  if(method == "MASTcpm"){
    message("Performing DEA using MAST with CPM data")
    timing <- system.time({
      names(group) <- colnames(data)
      dge <- edgeR::DGEList(counts = data)
      dge <- edgeR::calcNormFactors(dge)
      cpms <- edgeR::cpm(dge)
      sca <- MAST::FromMatrix(exprsArray = log2(cpms + 1),
                              cData = data.frame(wellKey = names(group),
                                                 group = group))
      zlmdata <- MAST::zlm(~ group, sca)
      mast <- MAST::lrTest(zlmdata, "group")
      df <- data.frame(PValue = mast[, "hurdle", "Pr(>Chisq)"],
                       FDR = stats::p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH"),
                       row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
    })
    return(df)
  }

  if(method == "MASTtpm"){
    message("Performing DEA using MAST with TPM data")
    tpm <- simutils::normalization_simutils(data = data,
                                            norm_method = "TPM",
                                            species = species)
    timing <- system.time({
      names(group) <- colnames(tpm)
      sca <- MAST::FromMatrix(exprsArray = log2(tpm + 1),
                              cData = data.frame(wellKey = names(group),
                                                 group = group))
      zlmdata <- MAST::zlm(~ group, sca)
      mast <- MAST::lrTest(zlmdata, "group")
      df <- data.frame(PValue = mast[, "hurdle", "Pr(>Chisq)"],
                       FDR = stats::p.adjust(mast[, "hurdle", "Pr(>Chisq)"], method = "BH"),
                       row.names = names(mast[, "hurdle", "Pr(>Chisq)"]))
    })
    return(df)
  }

  if(method == "limmatrend"){
    message("Performing DEA using limma-trend")
    timing <- system.time({
      dge <- edgeR::DGEList(data, group = group)
      dge <- edgeR::calcNormFactors(dge)
      design <- stats::model.matrix(~ group)
      y <- methods::new("EList")
      y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
      fit <- limma::lmFit(y, design = design)
      fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
      tt <- limma::topTable(fit, n = Inf, adjust.method = "BH")
    })
    return(tt$table)
  }

  if(method == "limmavoom"){
    message("Performing DEA using limma-voom")
    timing <- system.time({
      dge <- edgeR::DGEList(data, group = group)
      dge <- edgeR::calcNormFactors(dge)
      design <- model.matrix(~ group)
      vm <- limma::voom(dge, design = design, plot = FALSE)
      fit <- limma::lmFit(vm, design = design)
      fit <- limma::eBayes(fit)
      tt <- limma::topTable(fit, n = Inf, adjust.method = "BH")
    })
    return(tt$table)
  }

  if(method == "ttest"){
    message("Performing DEA using t-test")
    tpm <- simutils::normalization_simutils(data = data,
                                            norm_method = "TPM",
                                            species = species)
    timing <- system.time({
      tmm <- edgeR::calcNormFactors(tpm)
      tpmtmm <- edgeR::cpm(tpm, lib.size = tmm * colSums(tpm))
      logtpm <- log2(tpmtmm + 1)
      idx <- seq_len(nrow(logtpm))
      names(idx) <- rownames(logtpm)
      ttest_p <- sapply(idx, function(i) {
        stats::t.test(logtpm[i, ] ~ group)$p.value
      })
    })
    df <- data.frame(PValue = ttest_p,
                     FDR = stats::p.adjust(ttest_p, method = "BH"),
                     row.names = names(ttest_p))
    return(df)
  }

  if(method == "wilcox"){
    message("Performing DEA using wilcox-test")
    tpm <- simutils::normalization_simutils(data = data,
                                            norm_method = "TPM",
                                            species = species)
    timing <- system.time({
      tmm <- edgeR::calcNormFactors(tpm)
      tpmtmm <- edgeR::cpm(tpm, lib.size = tmm * colSums(tpm))
      idx <- 1:nrow(tpmtmm)
      names(idx) <- rownames(tpmtmm)
      wilcox_p <- sapply(idx, function(i) {
        stats::wilcox.test(tpmtmm[i, ] ~ group)$p.value
      })
    })
    df <- data.frame(PValue = wilcox_p,
                     FDR = stats::p.adjust(wilcox_p, method = "BH"),
                     row.names = names(wilcox_p))
    return(df)
  }
  print(timing)
}
