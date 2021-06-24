#' Estimate the gene dispersion
#'
#' @description Estimate the gene dispersion of each transcripts by using EdgeR estimateDisp function, then scale it to match OUTRIDER gene dispersion estimation.
#'
#' @param SequencingTable a read counts table with the transcripts in row and the samples in column.
#'
#' @return Vector of each transcripts gene dispersion.
#' @importFrom edgeR DGEList
#' @importFrom edgeR estimateDisp
#' @examples
#' SequencingTable <- ExampleAbeilleDataSet
#' genedisp <- GeneDisp(SequencingTable)
#' print(head(genedisp))
#' @export
GeneDisp <- function(SequencingTable) {
  y <- DGEList(SequencingTable)
  y <- estimateDisp(y)
  GeneDisp <- 1/(y$tagwise.dispersion)*10
  return(GeneDisp)
}
