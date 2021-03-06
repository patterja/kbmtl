% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_classification_celllines.R
\name{test_classification_celllines}
\alias{test_classification_celllines}
\title{DrugPrediction TESTING
version 1}
\usage{
test_classification_celllines(combined_matrix, train_idx, test_idx,
  cell_line_response = "cell_line_response_threshold_0.50_large_and_small_screen.RData",
  targetid = "target_id.txt", genelist = "Trusight_genes.csv", state)
}
\arguments{
\item{combined_matrix}{(string):path to FINAL_norm.txt run through normalizeRUV in SMMARTFunctions.R}

\item{train_idx}{(string):cell_line_response_threshold_0.50_large_and_small_screen.RData}

\item{test_idx}{(string):cell_line_response_threshold_0.50_large_and_small_screen.RData}

\item{cell_line_response}{(string): cell_line_response_threshold_0.50_large_and_small_screen.RData. Cell line are rows and drugs are columns.}

\item{targetid}{(string): path to target_id.txt. Converts to HUGO identifiers and filters to protein coding genes.}

\item{genelist}{(string): path to gene list Trusight.csv. Used to filter genes associated with cancer. Make sure gene identifiers match.}

\item{state}{(string): path to "smmart_trained_machine_learning_model.RData" from drugPredTrain}
}
\value{
Y_predicted(data.frame): prediction
}
\description{
DrugPrediction TESTING
version 1
}
\author{
Janice Patterson
}
