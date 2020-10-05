#' DrugPrediction TESTING
#
#' version 1
#' @param combined_matrix (string):path to FINAL_norm.txt run through normalizeRUV in SMMARTFunctions.R
#' @param train_idx (string):cell_line_response_threshold_0.50_large_and_small_screen.RData
#' @param test_idx (string):cell_line_response_threshold_0.50_large_and_small_screen.RData
#' @param cell_line_response (string): cell_line_response_threshold_0.50_large_and_small_screen.RData. Cell line are rows and drugs are columns.
#' @param targetid (string): path to target_id.txt. Converts to HUGO identifiers and filters to protein coding genes.  
#' @param genelist (string): path to gene list Trusight.csv. Used to filter genes associated with cancer. Make sure gene identifiers match. 
#' @param state (string): path to "smmart_trained_machine_learning_model.RData" from drugPredTrain
#' @return Y_predicted(data.frame): prediction
#' 
#' @export
#' @author Janice Patterson


test_classification_celllines <- function(combined_matrix,
                                          train_idx,
                                          test_idx,
                                          cell_line_response = "cell_line_response_threshold_0.50_large_and_small_screen.RData",
                                          targetid = "target_id.txt",
                                          genelist = "Trusight_genes.csv",
                                          state)
{
  version="0.1"
  
  if (is.null(combined_matrix)|is.null(train_idx)|is.null(cell_line_response)){print(parser$print_help())}
  if (is.null(combined_matrix)|is.null(train_idx)|is.null(cell_line_response)){print(parser$print_help())}
  
  
  #~ Ytrain ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  load(cell_line_response)
  Y_train <- cell_line_response
  
  #~ Xtrain ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  expression <- read.table(combined_matrix, row.names = 1, check.names = FALSE, sep = "\t", header=TRUE)
  #transpose matrix for matrix math rows are samples, genes are columns
  train_expression <- t(expression[,train_idx])
  
  #~ Xtest ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #Test Data
  test_expression <- t(expression[,test_idx])
  ## remove UHRs if they exist
  test_expression <- test_expression[grepl(rownames(test_expression), pattern = "UHR|horizon") == FALSE,]
  
  #~ YTRAIN FILTERING  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## filter sample: check to make sure samples of expression are included in samples of response
  common_cell_lines <- intersect(rownames(Y_train), rownames(train_expression))
  
  if (length(common_cell_lines)>0){
    train_expression = train_expression[common_cell_lines,]
    Y_train <- Y_train[common_cell_lines,]
  } else {
    stop("no samples in common with samples in response matrix and samples in expression matrix")
  }
  
  ## filter frequency: filtering for frequencies with a population of resistant and sensitive
  frequencies <- colSums(Y_train == 1, na.rm = TRUE) / colSums(is.na(Y_train) == FALSE)
  
  Y_train <- Y_train[,which(frequencies > 0.05 & frequencies < 0.95)]
  
  #~ XTRAIN FILTERING  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ## filter protein coding: filter for protein coding genes
  if (is.null(targetid)){
    print("Not filtering for protein coding genes.")
    ttrain_expression = t(train_expression)
    ttest_expression = t(test_expression)
  } else{
    ttrain_expression = convert2hugo(t(train_expression), targetid_file = targetid)
    ttest_expression = convert2hugo(t(test_expression), targetid_file = targetid)
    print("Filtering for protein coding genes and converting to HUGO names")}
  
  
  # FILTER trusight genes: filter gene lis of known cancer types
  if (is.null(genelist)){
    print("Not filtering by gene list. Specify list of cancer related genes to train on.")
    ttrain_expression = ttrain_expression
    ttest_expression = ttest_expression
  } else {
    trusight_genes <- read.csv(genelist, header = FALSE, stringsAsFactors = FALSE)[,1]
    ltru = intersect(trusight_genes, rownames(ttrain_expression))
    ttrain_expression <- ttrain_expression[ltru,]
    ttest_expression <- ttest_expression[ltru,]
    }
  
  
  ltrain_expression = log2(t(ttrain_expression + 1))
  ltest_expression = log2(t(ttest_expression +1))

  X_train <- scale(ltrain_expression)
  X_test <- scale(ltest_expression)
  
  valid_genes <- which(colSums(is.na(X_train)) == 0)
  X_train <- X_train[, valid_genes]
  X_test <- X_test[, valid_genes]
  X_test[is.na(X_test)] <- 0
  

  K_train <- X_train %*% t(X_train)
  K_test <- X_train %*% t(X_test)
  normalizer <- max(abs(K_train))
  K_train <- K_train / normalizer
  K_test <- K_test / normalizer

  #training
  load(state)
  prediction <- kbmtl_semisupervised_classification_variational_test(K_test,state)

  Y_predicted <- as.data.frame(prediction$P)
  colnames(Y_predicted) <- colnames(Y_train)
  return(Y_predicted)
}
