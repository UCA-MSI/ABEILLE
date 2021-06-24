#' Estimate the size factor
#'
#' @description Estimate the size factor of each samples by using DESeq2 estimateSizeFactors function.
#'
#' @param SequencingTable a read counts table with the transcripts in row and the samples in column.
#' @param annotation_table a table of containing characteristics on each samples.
#'
#' @return Vector where each elements is the size factor of a sample.
#' @importFrom DESeq2 DESeqDataSetFromMatrix
#' @importFrom DESeq2 estimateSizeFactors
#'
#' @examples
#' anno_table <- data.frame(sampleID = colnames(ExampleAbeilleDataSet),
#' Sex = sample(c("Male", "Female"), 128, replace = TRUE),
#' Age = sample(1:89, 128, replace = TRUE), Cohort = rep("GTEx database",128))
#' SequencingTable <- ExampleAbeilleDataSet
#' sizefactor <- SizeFactors(ExampleAbeilleDataSet, anno_table)
#' print(head(sizefactor))
#' @export
SizeFactors <- function(SequencingTable, annotation_table) {
  dds <- DESeqDataSetFromMatrix(countData = SequencingTable, colData = annotation_table, design = ~1)
  dds <- estimateSizeFactors(dds)
  return(dds@colData@listData[["sizeFactor"]])
}
