#' Perform Differential Expression Analysis
#' This function is used to perform DEA for single-cell RNA-seq data by QLF model in edgeR package including the cellular detection rate.
#' @param data A matrix of gene expression profile.
#' @param group Group(or cluster) assignment of every cells in columns of matrix.
#' @param method Which DE method will be used?
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest topTags cpm
#' @importFrom stats model.matrix
#' @return A dataframe of DEA result.
#' @export
#' @references
#' Soneson C, Robinson M D. Bias, robustness and scalability in single-cell differential expression analysis. Nature methods, 2018, 15(4): 255-261. <https://doi.org/10.1038/nmeth.4612>
perform_DEA <- function(
  data,
  group,
  method
){
  ## Filter
  filter_index <- rowSums(data) > 0
  message(paste0(nrow(data) - sum(filter_index), " genes are removed when filtering"))
  data <- data[rowSums(data) > 0, ]

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
  }
  if(method == "limmatrend"){
    message("Performing DEA using limma-trend")
    timing <- system.time({
      dge <- edgeR::DGEList(data, group = group)
      dge <- edgeR::calcNormFactors(dge)
      design <- stats::model.matrix(~ group)
      y <- new("EList")
      y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
      fit <- limma::lmFit(y, design = design)
      fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
      tt <- edgeR::topTable(fit, n = Inf, adjust.method = "BH")
    })
  }

  if(method == "limmavoom"){
    timing <- system.time({
      dge <- edgeR::DGEList(data, group = group)
      dge <- edgeR::calcNormFactors(dge)
      design <- model.matrix(~ group)
      vm <- limma::voom(dge, design = design, plot = FALSE)
      fit <- limma::lmFit(vm, design = design)
      fit <- limma::eBayes(fit)
      tt <- edgeR::topTable(fit, n = Inf, adjust.method = "BH")
    })
  }

  print(timing)
  result <- tt$table
  return(result)
}



# cd <- clean.counts(data, min.lib.size=100, min.reads = 3, min.detected = 3)
# group <- group[names(group) %in% colnames(cd)]
#
#
# o.ifm <- scde.error.models(counts = cd,
#                            groups = group,
#                            n.cores = 8,
#                            threshold.segmentation = TRUE,
#                            save.crossfit.plots = FALSE,
#                            save.model.plots = FALSE,
#                            verbose = 1)
# head(o.ifm)
#
# valid.cells <- o.ifm$corr.a > 0
# table(valid.cells)
#
# o.ifm <- o.ifm[valid.cells, ]
#
# # estimate gene expression prior
# o.prior <- scde.expression.prior(models = o.ifm,
#                                  counts = cd,
#                                  length.out = 400,
#                                  show.plot = FALSE)
# groups <- factor(group[match(rownames(o.ifm), names(group))], levels = c("PBS", "LPS"))
# names(groups) <- rownames(o.ifm)
# ediff <- scde.expression.difference(o.ifm,
#                                     cd,
#                                     o.prior,
#                                     groups  =  groups,
#                                     n.randomizations  =  100,
#                                     n.cores  =  6,
#                                     verbose  =  1)
# head(ediff[order(ediff$Z, decreasing  =  TRUE), ])
# ediff$pvalue <- 2*pnorm(-abs(ediff$Z))
#
# batch <- as.factor(ifelse(rbinom(nrow(o.ifm), 1, 0.5) == 1, "batch1", "batch2"))
# table(groups, batch)
#
# ediff.batch <- scde.expression.difference(o.ifm,
#                                           cd,
#                                           o.prior,
#                                           groups = groups,
#                                           batch = batch,
#                                           n.randomizations = 100,
#                                           n.cores = 4,
#                                           return.posteriors = TRUE,
#                                           verbose = 1)
#
# library(ROTS)
# data(upsSpikeIn)
# input = upsSpikeIn
# groups = c(rep(0,3), rep(1,3))
# results = ROTS(data = input, groups = groups , B = 100 , K = 500 , seed = 1234)
#
# t1 <- Sys.time()
# result <- ROTS(data = data, groups = group, seed = 111)
# t2 <- Sys.time()
# difftime(t2, t1)
# de_table <- data.frame("pvalue" = result$pvalue,
#                        "FDR" = result$pvalue,
#                        "logfc" = result$logfc)
# de_gene_rots <- de_table %>% filter(FDR < 0.05)
