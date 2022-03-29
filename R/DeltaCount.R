#' Compute the Delta Count
#'
#' @description Metric used to find AGE, this definition is inspired from the log2 fold-change is not the classic one.
#'
#' @param SequencingTable a read counts table with the transcripts in row and the samples in column.
#' @param ReconstructedTable the reconstructed table generate by an autoencoder.
#'
#' @return All Delta Count in a matrix. The computation is the following one:
#' \ifelse{html}{
#'     \out{Δ<sub>ij</sub> = log2(k<sub>ij</sub> / μ<sub>i</sub><sup>k&#770)}}{
#'     \deqn{Δ_{ij} = log2(k_{ij} / μ&hat_{ij}^{k&hat})}}.
#' More detail in the ABEILLE paper.
#' @examples
#' SequencingTable <- ExampleAbeilleDataSet
#' ReconstructedTable <- ExampleAbeilleReconstructed
#' delta_count <- DeltaCount(SequencingTable, ReconstructedTable)
#' print(head(delta_count))
#' @export
DeltaCount <- function(SequencingTable, ReconstructedTable) {
  df_to_fill <- data.frame(matrix(nrow=nrow(SequencingTable),ncol=ncol(SequencingTable)))
  list_col <- 1:ncol(SequencingTable)
  df_to_fill[,1] <- rowMeans(SequencingTable[,list_col[-1]])
  for (i in 2:length(list_col)) {
    df_to_fill[,i] <- rowMeans(SequencingTable[,list_col[-i]])
  }
  rownames(df_to_fill) <- rownames(SequencingTable)
  colnames(df_to_fill) <- colnames(SequencingTable)
  delta_count <- log2(ReconstructedTable/df_to_fill)
  delta_count <- as.matrix(delta_count)
  delta_count[!is.finite(delta_count)] <- 0
  return(delta_count)
}
