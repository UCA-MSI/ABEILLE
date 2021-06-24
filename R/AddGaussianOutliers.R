#' Add outliers using gaussian noise
#'
#' @description Create a copy of a sequencing table with outliers using gaussian noise.
#'
#' @param SequencingTable a read counts table with the transcripts in row and the samples in column.
#' @param drawoption the way to pick outliers, if the value is smaller than 1 it is the probability of observation gets noisy. If it is higher than one it is the number of observations to change.
#' @param feedback logical, add a data frame containing informations on the outliers.
#' @param setseed set the seed.
#'
#' @return A sequencing table with outliers. Outliers are made using the formula:
#' \ifelse{html}{
#'     \out{k<sub>ij</sub><sup>&#119977</sup> = k<sub>ij</sub> + &#119977(0,σ<sub>i</sub>)}}{
#'     \deqn{k_{ij}^{N}=k_{ij}+N(0,σ_{j}^{N})}}.
#' More detail in the ABEILLE paper.
#' @importFrom  stats rbinom
#' @importFrom stats rnorm
#' @examples
#' SequencingTable <- ExampleAbeilleDataSet
#' TableWithOutliers <- AddGaussianOutliers(SequencingTable,
#' drawoption = 10^-4, feedback = TRUE, setseed = 5)
#' print(head(TableWithOutliers$table))
#' print(head(TableWithOutliers$feedback))
#' @export
AddGaussianOutliers <- function(SequencingTable, drawoption = 10^-5, feedback=FALSE, setseed=FALSE){
  if (class(SequencingTable) == "data.frame"){SequencingTable <- as.matrix(SequencingTable)}
  if (setseed != FALSE){
    set.seed(setseed)
  }
  if (drawoption < 1){
    nbdraw <- rbinom(1, dim(SequencingTable)[1] * dim(SequencingTable)[2], drawoption)
  } else{
    nbdraw <- drawoption
  }
  random_row <- sample(1:dim(SequencingTable)[1], nbdraw, replace = TRUE)
  random_col <- sample(1:dim(SequencingTable)[2], nbdraw, replace = TRUE)
  TableWithOutliers <- SequencingTable
  original_value <- c()
  transformed_value <- c()
  for (i in 1:length(random_row)){
    original_value <- c(original_value,SequencingTable[random_row[i],random_col[i]])
    TableWithOutliers[random_row[i],random_col[i]] <- TableWithOutliers[random_row[i],random_col[i]] +  rnorm(1,0,sd(TableWithOutliers[random_row[i],]))
    transformed_value <- c(transformed_value,TableWithOutliers[random_row[i],random_col[i]])
  }
  if (feedback == FALSE){
    return(TableWithOutliers)
  } else {
    feedbackdata <- data.frame("Sample" = colnames(SequencingTable)[random_col],"Transcript" = rownames(SequencingTable)[random_row],"Original" = original_value, "Transformed" = transformed_value)
    listreturn <- list('table' = TableWithOutliers, 'feedback' = feedbackdata)
    return(listreturn)
  }
}
