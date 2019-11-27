#' Train KBMTL classifer on Breast Cancer celllines
#' 
#' @param combined_matrix (string): Tab separated table containing training and test data. Genes are in the rows and samples are in the columns. 
#' @param train_idx (numeric): Numeric comma separated list of index of the columns that contain the training data from the table in path_to_combined. Index 0 is the gene names. Ranges are colon spearated. eg. 1,2,4:10
#' @param cell_line_response (string): cell_line_response_threshold_0.50_large_and_small_screen.RData. Cell line are rows and drugs are columns.
#' @param targetid (string): path to target_id.txt. Converts to HUGO identifiers and filters to protein coding genes.  
#' @param genelist (string): path to gene list Trusight.csv. Used to filter genes associated with cancer. Make sure gene identifiers match. 
#' @return smmart_trained_machine_learning_model.RData
#' 
#' @export
#' @author Janice Patterson
#' 



train_classification_celllines <- function(combined_matrix, 
                                           train_idx, 
                                           cell_line_response = "cell_line_response_threshold_0.50_large_and_small_screen.RData",
                                           targetid = "target_id.txt",
                                           genelist = "Trusight_genes.csv"){

version="0.1"

if (is.null(combined_matrix)|is.null(train_idx)|is.null(cell_line_response)){print(parser$print_help())}

  
#~ Ytrain ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load(cell_line_response)
Y_train <- cell_line_response

#~ Xtrain ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
expression <- read.table(combined_matrix, row.names = 1, check.names = FALSE, sep = "\t", header=TRUE)
#transpose matrix for matrix math rows are samples, genes are columns
train_expression <- t(expression[,train_idx])


# Ytrain FILTERING  ~~~~~~~~~~~~~~~~~~~~~~~~
# FILTER sample: check to make sure samples of expression are included in samples of response
common_cell_lines <- intersect(rownames(Y_train), rownames(train_expression))

if (length(common_cell_lines)>0){
  train_expression = train_expression[common_cell_lines,]
  Y_train <- Y_train[common_cell_lines,]
  } else {
  stop("no samples in common with samples in response matrix and samples in expression matrix")
}

# FILTER frequency: filtering for frequencies with a population of resistant and sensitive
frequencies <- colSums(Y_train == 1, na.rm = TRUE) / colSums(is.na(Y_train) == FALSE)

Y_train <- Y_train[,which(frequencies > 0.05 & frequencies < 0.95)]

# Xtrain FILTERING  ~~~~~~~~~~~~~~~~~~~~~~~~
# FILTER protein coding: filter for protein coding genes
if (is.null(targetid)){
  print("Not filtering for protein coding genes.")
  ttrain_expression = t(train_expression)
  } else{
  ttrain_expression = convert2hugo(t(train_expression), targetid_file = targetid)
  print("Filtering for protein coding genes and converting to HUGO names")}


# FILTER trusight genes: filter gene lis of known cancer types
if (is.null(genelist)){
  print("Not filtering by gene list. Specify list of cancer related genes to train on.")
  ttrain_expression = ttrain_expression
  } else {
  trusight_genes <- read.csv(genelist, header = FALSE, stringsAsFactors = FALSE)[,1]
  ltru = intersect(trusight_genes, rownames(ttrain_expression))
  ttrain_expression <- ttrain_expression[ltru,]
}


ltrain_expression = log2(t(ttrain_expression + 1))

X_train <- scale(ltrain_expression)

#FILTER remove genes that had no variation 
valid_genes <- which(colSums(is.na(X_train)) == 0)
X_train <- X_train[, valid_genes]
  

  
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
