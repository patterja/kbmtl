#' Train KBMTL classification model generic
#' 
#' @param Xtrain (matrix): training data. Featrues are columns and samples are rows 
#' @param Ytrain (matrix): response 
#' @return smmart_trained_machine_learning_model.RData
#' @export
#' @example 
#' 




train_classification <- function(Xtrain, Ytrain) {
  #parse optional arguments 
  
  
  K_train <- X_train %*% t(X_train)
  normalizer <- max(abs(K_train))
  K_train <- K_train / normalizer
  
  parameters <- list()
  parameters$alpha_lambda <- 1
  parameters$beta_lambda <- 1
  parameters$iteration <- 200
  parameters$margin <- 1
  parameters$R <- 20
  parameters$seed <- 1606
  parameters$sigma_h <- 0.1
  parameters$sigma_w <- 1.0
  
  state <- kbmtl_semisupervised_classification_variational_train(K_train, Y_train, parameters)
  
  save(state, file =  "smmart_trained_machine_learning_model.RData")
  print(paste0("kbmtl trained model saved ", getwd(), "/smmart_trained_machine_learning_model.RData"))
  return(state)
}
