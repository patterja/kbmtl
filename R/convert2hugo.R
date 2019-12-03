#' convert ENSGids to HUGO gene names 
#' 
#' RNA matrix of ENSG identifiers converted and aggregated (summed) to HUGO
#' identifiers using targetid file from GENCODE v24. 
#' @param mat (matrix) columns are sample names, rows are ENSG IDs
#' @param targetid_file (char string) path to target_id.txt
#' @param protein_coding (logical): filter for protein coding only and then aggregate to HUGO
#' @return aggregated matrix with HUGO IDs
#' @export
#' @author Janice Patterson

convert2hugo <- function(mat, targetid_file="target_id.txt", protein_coding=TRUE){
  targetid = read.csv(targetid_file, sep="\t", stringsAsFactors = F)
  
  mat=as.matrix(mat)
  #filter by for protein coding
  if (protein_coding==TRUE){
    targetid <- targetid[which(targetid$X.7 == "protein_coding"),]}
  #remove redundant combinations of HUGO and ENSG
  targetid <- unique(targetid[, c("X.1", "X.5")])
  rownames(targetid) <- targetid$X.1
  #filter the mat for the targets
  mat = as.matrix(mat)
  mat = mat[intersect(targetid[,"X.1"], rownames(mat)),]
  #convert to HUGO
  rownames(mat) = targetid$X.5[match(rownames(mat), targetid$X.1)]
  #redundant HUGO names are summed
  mat_gene = aggregate(mat,by=list(rownames(mat)), FUN=sum)
  row.names(mat_gene) <- mat_gene[[1]]
  mat_gene[[1]]=NULL

  return(mat_gene)
}
