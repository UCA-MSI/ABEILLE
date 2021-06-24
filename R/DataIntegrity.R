#' Check the integrity of a DataSet
#'
#' @description Verifies that the data is usable and print the potential problems.
#'
#' @param SequencingTable a read counts table with the transcripts in row and the samples in column.
#'
#' @return Printing message about potential dataset issues.
#' @importFrom crayon black
#' @importFrom crayon green
#' @examples
#' SequencingTable <- ExampleAbeilleDataSet
#' DataIntegrity(SequencingTable)
#' @export
DataIntegrity <- function(SequencingTable) {
  if(length(dim(SequencingTable)) != 2)
    stop("Data dimension is incorrect, expected 2D array")
  if(dim(SequencingTable)[1] < dim(SequencingTable)[2])
    cat(black("Check that the transcripts  are in rows and the samples in columns\n"))
  if(any(is.na(SequencingTable)))
    stop("The data contains NA(s) value(s)")
  if(any(is.infinite(as.matrix(SequencingTable))))
    stop("The data contains infinite(s) value(s)")
  if(any(is.factor(SequencingTable)))
    stop("The data contains factor(s) value(s)")
  if(any(!apply(SequencingTable,2,is.numeric)))
    stop("Some value(s) is/are not numeric")
  cat(green("Data integrity is validated"))
}
