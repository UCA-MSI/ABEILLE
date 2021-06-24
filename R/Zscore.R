#' Compute the Z-score
#'
#' @description Metric used to find outliers, this definition of the Z-score is not the classic one.
#'
#' @param SequencingTable a read counts table with the transcripts in row and the samples in column.
#' @param ReconstructedTable the reconstructed table generate by an autoencoder.
#'
#' @return All log2 fold-change in a matrix. The computation is the following one:
#' \ifelse{html}{
#'     \out{L<sub>ij</sub> = log2((k<sub>ij</sub> + 1) / (k&#770<sub>ij</sub> + 1)) and Z<sub>ij</sub> = (L<sub>ij</sub> - μ<sub>i</sub><sup>L</sup>)/σ<sub>i</sub><sup>L</sup>}}{
#'     \deqn{L_{ij} = log2((k_{ij}+ 1) / (k&hat_{ij}+ 1)) and Z_{ij}> = (L_{ij} - μ_{i}^{L})/σ_{i}^{L}}}.
#' More detail in the ABEILLE paper.
#' @importFrom matrixStats rowSds
#' @export
Zscore <- function(SequencingTable, ReconstructedTable){
  if (class(SequencingTable)[1] == "data.frame"){SequencingTable <- as.matrix(SequencingTable)}
  if (class(ReconstructedTable)[1] == "data.frame"){ReconstructedTable <- as.matrix(ReconstructedTable)}
  Lij <- log2((SequencingTable+1)/(ReconstructedTable+1))
  Lij <- (Lij - rowMeans(Lij)) / rowSds(Lij)
  Lij[!is.finite(Lij)] <- 0
  return(Lij)
}
