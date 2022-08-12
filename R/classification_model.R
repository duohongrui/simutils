# pre_data <- readRDS("/Users/duohongrui/Desktop/preprocessed_data/data1_GSE54006.rds")
# data <- pre_data$data
# data_info <- pre_data$data_info
# group <- data_info$treatment
#
# result <- simutils::perform_DEA(data = data, group = group, method = "edgeRQLF")
# de_genes <- rownames(result)[result$FDR < 0.05]


#' Establish A Model And Predict Using DEGs
#'
#' @param data A matrix with cells on columns and genes on rows.
#' @param group Group(or cluster) assignment of every cells in columns of matrix.
#' @param de_genes A character vector of DEGs.
#' @param method The method to establish the model. SVM, Decision tree.
#' @importFrom stats predict
#' @export
model_predict <- function(
  data,
  group,
  de_genes,
  method
){
  ## filter genes with constant value across cells
  gene_var <- apply(data, 1, BiocGenerics::var)
  data <- data[gene_var != 0, ]
  ## scale
  data <- scale(data, center = FALSE)
  ## split data into train and test subsets
  set.seed(1111)
  train_index <- sample(1:ncol(data), ncol(data) * 0.8, replace = FALSE)
  data <- as.data.frame(t(data))
  train_data <- data[train_index, intersect(de_genes, names(gene_var)[gene_var != 0])]
  test_data <- data[-train_index, intersect(de_genes, names(gene_var)[gene_var != 0])]
  ## group information
  group <- as.numeric(as.factor(group))
  train_group <- group[train_index]
  test_group <- group[-train_index]

  if(method == "SVM"){
    ## SVM
    if(!requireNamespace("e1071", quietly = TRUE)){
      stop("Package \"e1071\" must be installed by \"install.packages('e1071')\" command.")
    }
    svm_classifier <- e1071::svm(x = train_data,
                                 y = train_group,
                                 cross = 10,
                                 kernel = 'radial',
                                 prob = TRUE,
                                 scale = FALSE)
    svm_predict <- stats::predict(svm_classifier,
                                  as.matrix(test_data),
                                  prob = TRUE)
    return(svm_predict)
  }

  if(method == "Decision tree"){
    ## Decision tree
    if(!requireNamespace("rpart", quietly = TRUE)){
      stop("Package \"rpart\" must be installed by \"install.packages('rpart')\" command.")
    }
    train_data$group <- train_group
    tree_model <- rpart::rpart(group ~ .,
                               data = train_data,
                               method = "class")
    tree_predict <- stats::predict(tree_model,
                                   as.data.frame(test_data),
                                   prob = TRUE)
    return(tree_predict)
  }
}



# par(pty = "s")
# print(pROC::roc(test_group,
#                 tree_predict[,2],
#                 plot = TRUE,
#                 legacy.axes = TRUE,
#                 percent = TRUE,
#                 col = "#377eb8",
#                 print.auc = TRUE,
#                 xlab = "FPR",
#                 ylab = "TPR",
#                 lwd = 3))
# roc_obj <- pROC::roc(test_group, svm_predict)
