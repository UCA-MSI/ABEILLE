#' Remove zero counts
#'
#' @description Function to remove all transcripts with 0 count.
#'
#' @param SequencingTable a read counts table with the transcripts in row and the samples in column.
#'
#' @return The table with the 0 count row removed and a print to indicate the number of transcript deleted.
#' @examples
#' SequencingTable <- ExampleAbeilleDataSet
#' WithoutZeroCounts <- RemoveZeroCounts(SequencingTable)
#' print(head(WithoutZeroCounts))
#' @export
RemoveZeroCounts <- function(SequencingTable) {
  cat(nrow(SequencingTable) - sum((rowSums(SequencingTable) != 0)),
      "transcripts removed","over",
      nrow(SequencingTable),
      "(which corresponds to",
      round(sum((rowSums(SequencingTable) == 0))/nrow(SequencingTable),4)*100,
      "% of the dataset)")
  SequencingTable <- SequencingTable[(rowSums(SequencingTable) != 0),]
  return(SequencingTable)
}
