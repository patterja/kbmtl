#' Training functions for classification on cell lines with prediction and expression
#'
#' @param X_train (matrix) matrix of expression
#' @param Y_train (matrix) matrix of response binary (0,1)
#' @return state
#' @export
#' @author Janice Patterson
#' @examples state = cellline_train(X_train, Y_train)

cellline_train <- function(X_train, Y_train) {
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

  state <- kbmtl_semisupervised_classification_variational_train(K_train,
Y_train, parameters)
  return(state)
}

