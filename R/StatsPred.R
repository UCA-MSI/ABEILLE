#' Display statistics on the error between of the input and the output of an autoencoder
#'
#' @description Print in the console different statistics of the absolute error between the two datasets provided.
#'
#' @param SequencingTable a read counts table with the transcripts in row and the samples in column.
#' @param ReconstructedTable the reconstructed table generate by an autoencoder.
#'
#' @return Print in the console these statistics.
#' @importFrom crayon bold
#' @importFrom stats median
#' @importFrom stats quantile
#' @importFrom stats sd
#' @importFrom e1071 skewness
#'
#' @examples
#' SequencingTable <- ExampleAbeilleDataSet
#' ReconstructedTable <- ExampleAbeilleReconstructed
#' StatsPred(SequencingTable,ReconstructedTable)
#' @export
StatsPred <- function(SequencingTable, ReconstructedTable) {
  if (class(SequencingTable) == "data.frame"){SequencingTable <- as.matrix(SequencingTable)}
  if (class(ReconstructedTable) == "data.frame"){ReconstructedTable <- as.matrix(ReconstructedTable)}
  list_dist <- abs(SequencingTable - ReconstructedTable)
  return(
    cat(
      bold("Number of observations :"), ncol(ReconstructedTable) * nrow(ReconstructedTable),
      bold("\nMean :"), mean(list_dist),
      bold("\nMedian :"), median(list_dist),
      bold("\nMin :"), min(list_dist),
      bold("\nQuantile 10% :"), quantile(list_dist,0.1),
      bold("\nQuantile 25% :"), quantile(list_dist,0.25),
      bold("\nQuantile 75% :"), quantile(list_dist,0.75),
      bold("\nQuantile 90% :"), quantile(list_dist,0.90),
      bold("\nMax :"), max(list_dist),
      bold("\nStandard deviation :"), sd(list_dist),
      bold("\nSkewness :"), skewness(list_dist),
      bold("\nEmpty samples:"), sum(colSums(ReconstructedTable) == 0),bold("on"),dim(ReconstructedTable)[2],
      bold("\nEmpty transcripts:"), sum(rowSums(ReconstructedTable) == 0),bold("on"),dim(ReconstructedTable)[1],
      bold("\nNumber of zeros:"), sum(ReconstructedTable == 0), bold("on"), sum(SequencingTable == 0)
    )
  )
}
