#' test generic
#' 
#' #' Testing function for classification on cell lines with prediction and expression
#' @param X_train (matrix) matrix of expression
#' @param X_test (matrix) matrix
#' @param Y_train (matrix) matrix of response binary (0,1)
#' @param state (vector)
#' @return prediction
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
