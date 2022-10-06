# pre_data <- readRDS("/Users/duohongrui/Desktop/preprocessed_data/data1_GSE54006.rds")
# data <- pre_data$data
# data_info <- pre_data$data_info
# group <- data_info$cluster_info
#
# result <- simutils::perform_DEA(data = data, group = group, method = "edgeRQLF")
# de_genes <- rownames(result)[result$FDR < 0.05]
#
#
# estimate_result <- simmethods::Splat_estimation(ref_data = data, verbose = TRUE, seed = 111)
# prob.group = as.numeric(table(group)/ncol(data))
# de.prob = length(de_genes)/nrow(data)
# simulate_result <- simmethods::Splat_simulation(
#   estimate_result$estimate_result,
#   other_prior = list(prob.group = c(0.3, 0.4, 0.3),
#                      de.prob = 0.3),
#   return_format = "list",
#   verbose = TRUE,
#   seed = 111)
#
# data <- simulate_result$simulate_result$count_data
# col_data <- simulate_result$simulate_result$col_meta
# group <- col_data$group
# row_data <- simulate_result$simulate_result$row_meta
# de_genes <- rownames(row_data)[row_data$de_gene == "yes"]
# non_gene_genes <- rownames(row_data)[row_data$de_gene == "no"]
#
# data <- data[non_gene_genes, ]



#' Establish A Model And Predict Cell Identity Using DEGs
#'
#' @param data A matrix with cells on columns and genes on rows.
#' @param group Group(or cluster) assignment of every cells in columns of matrix.
#' @param de_genes A character vector of DEGs.
#' @param method The method to establish the model. SVM, Decision tree or RF (Random Forest).
#' @importFrom stats predict
#' @importFrom caret confusionMatrix
#' @importFrom pROC roc multiclass.roc
#' @export
model_predict <- function(
  data,
  group,
  de_genes,
  method
){
  ## filter genes with constant value across cells
  message("Preprocessing data...")
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
  group <- as.factor(group)
  train_group <- group[train_index]
  test_group <- group[-train_index]

  if(method == "SVM"){
    ## SVM
    if(!requireNamespace("e1071", quietly = TRUE)){
      message("e1071 is not installed on your device...")
      message("Installing e1071...")
      utils::install.packages("e1071")
    }
    message("Modeling by SVM...")
    svm_classifier <- e1071::svm(x = train_data,
                                 y = as.factor(train_group),
                                 cross = 10,
                                 probability = TRUE,
                                 kernel = 'radial',
                                 scale = FALSE)
    message("Predicting...")
    predict_class <- stats::predict(svm_classifier,
                                    as.matrix(test_data),
                                    prob = FALSE)
    conf_matrix <- caret::confusionMatrix(predict_class,
                                          test_group,
                                          mode = "everything")
    predict_prob <- stats::predict(svm_classifier,
                                   as.matrix(test_data),
                                   prob = TRUE)
    if(nlevels(group) == 2){
      roc <- pROC::roc(response = test_group,
                       predictor = attr(predict_prob, "probabilities"))
    }else{
      roc <- pROC::multiclass.roc(response = test_group,
                                  predictor = attr(predict_prob, "probabilities"))
    }
  }

  if(method == "Decision tree"){
    ## Decision tree
    if(!requireNamespace("rpart", quietly = TRUE)){
      message("rpart is not installed on your device...")
      message("Installing rpart...")
      utils::install.packages("rpart")
    }
    train_data$group <- train_group
    tree_model <- rpart::rpart(group ~ .,
                               data = train_data,
                               method = "class")
    predict_class <- stats::predict(tree_model,
                                    as.data.frame(test_data),
                                    type = "class")
    conf_matrix <- caret::confusionMatrix(predict_class,
                                          test_group,
                                          mode = "everything")
    predict_prob <- stats::predict(tree_model,
                                   as.data.frame(test_data),
                                   type = "prob")
    if(nlevels(group) == 2){
      roc <- pROC::roc(response = test_group,
                       predictor = predict_prob)
    }else{
      roc <- pROC::multiclass.roc(response = test_group,
                                  predictor = predict_prob)
    }
  }

  if(method == "RF"){
    if(!requireNamespace("randomForest", quietly = TRUE)){
      message("randomForest is not installed on your device...")
      message("Installing randomForest...")
      utils::install.packages("randomForest")
    }
    train_data$group <- train_group
    colnames(train_data) <- gsub(pattern = "-",
                                 replacement = "",
                                 colnames(train_data))
    colnames(test_data) <- gsub(pattern = "-",
                                 replacement = "",
                                 colnames(test_data))
    RF_model <- randomForest::randomForest(as.factor(group) ~ .,
                                           data = train_data)
    predict_class <- stats::predict(RF_model,
                                    test_data,
                                    type = "class")
    conf_matrix <- caret::confusionMatrix(predict_class,
                                          test_group,
                                          mode = "everything")
    predict_prob <- stats::predict(RF_model,
                                   test_data,
                                   type = "prob")
    if(nlevels(group) == 2){
      roc <- pROC::roc(response = test_group,
                       predictor = predict_prob)
    }else{
      roc <- pROC::multiclass.roc(response = test_group,
                                  predictor = predict_prob)
    }
  }
  return(dplyr::lst(conf_matrix, roc))
}


# par(pty = "s")
# print(pROC::multiclass.roc(result[["roc"]][["response"]],
#                 result[["roc"]][["predictor"]],
#                 plot = TRUE,
#                 legacy.axes = TRUE,
#                 percent = TRUE,
#                 col = "#377eb8",
#                 print.auc = TRUE,
#                 xlab = "FPR",
#                 ylab = "TPR",
#                 lwd = 3))
# roc_obj <- pROC::roc(test_group, RF_predict_prob[, 2])
