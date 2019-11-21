#' DrugPrediction TRAINING
#' 
#' @param path_to_bccl (string): Path to BCCL normCPMRUV_BCCL.tsv
#' @param path_to_test (string): Path to Test matrix 
#' @param path_to_combined (string): path to combined training and test matrix
#' @param train_idx (numeric): numeric list of index of training data sets in path_to_combined
#' @param cell_line_response (string): cell_line_response_threshold_0.50_large_and_small_screen.RData
#' @param path_to_targetid (string): path to target_id.txt. Used to filter protein_coding genes only. 
#' @param path_to_trusight (string): path to gene list Trusight.csv. Used to filter genes associated with cancer
#' @return smmart_trained_machine_learning_model.RData
#' @export
#'
#' @author Janice Patterson
#' @examples Y_predicted = cellline_test(Xtrain, X_test, Ytrain, state) 

cellline_test <- function(Xtrain, X_test, Ytrain, state)
{

  K_train <- X_train %*% t(X_train)
  K_test <- X_train %*% t(X_test)
  normalizer <- max(abs(K_train))
  K_train <- K_train / normalizer
  K_test <- K_test / normalizer

  #training
  prediction <- kbmtl_semisupervised_classification_variational_test(K_test,
state)

  Y_predicted <- as.data.frame(prediction$P)
  colnames(Y_predicted) <- colnames(Y_train)
  return(Y_predicted)
}
