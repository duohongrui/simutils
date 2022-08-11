# pre_data <- readRDS("/Users/duohongrui/Desktop/preprocessed_data/data1_GSE54006.rds")
# data <- pre_data$data
# data_info <- pre_data$data_info
# group <- data_info$treatment

# model_predict <- function(
#   data,
#   group,
#   de_genes,
#   svm = TRUE
# ){
#   ## filter genes with constant value across cells
#   gene_var <- apply(data, 1, BiocGenerics::var)
#   data <- data[gene_var != 0, ]
#   ## scale
#   data <- scale(data, center = FALSE)
#   ## split data into train and test subsets
#   set.seed(1111)
#   train_index <- sample(1:ncol(data), ncol(data) * 0.8, replace = FALSE)
#   data <- as.data.frame(t(data))
#   train_data <- data[train_index, intersect(de_genes, names(gene_var)[gene_var != 0])]
#   test_data <- data[-train_index, intersect(de_genes, names(gene_var)[gene_var != 0])]
#   ## group information
#   group <- as.numeric(as.factor(group))
#   train_group <- group[train_index]
#   test_group <- group[-train_index]
#
#   if(SVM){
#     ## SVM
#     svm_classifier <- e1071::svm(x = train_data,
#                                  y = train_group,
#                                  cross = 10,
#                                  kernel = 'radial',
#                                  prob = TRUE,
#                                  scale = FALSE)
#     svm_predict <- predict(svm_classifier,
#                            as.matrix(test_data),
#                            prob = TRUE)
#     return(svm_predict)
#   }
# }



# par(pty = "s")
# print(pROC::roc(test_group,
#                 svm_predict,
#                 plot = TRUE,
#                 legacy.axes = TRUE,
#                 percent = TRUE,
#                 col = "#377eb8",
#                 print.auc = TRUE,
#                 xlab = "FPR",
#                 ylab = "TPR",
#                 lwd = 3))
# roc_obj <- pROC::roc(test_group, svm_predict)
