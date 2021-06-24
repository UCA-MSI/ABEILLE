#' Add outliers using OUTRIDER method
#'
#' @description Create a copy of a sequencing table with outliers using OUTRIDER formula.
#'
#' @param SequencingTable a read counts table with the transcripts in row and the samples in column.
#' @param SizeFactor vector containing all the size factor. The size factor of each sample can be compute using the function SizeFactors, the function use the estimatation of DESeq2.
#' @param Zmean mean of the gaussian distribution used. Default log(3).
#' @param Zsd standard deviation of the gaussian distribution used. Default log(1.6).
#' @param drawoption the way to pick outliers, if the value is smaller than 1 it is the probability of observation gets noisy. If it is higher than one it is the number of observations to change.
#' @param feedback logical, add a data frame containing informations on the outliers.
#' @param setseed set the seed.
#'
#' @return A sequencing table with outliers. Outliers are drawn using the formula:
#' \ifelse{html}{
#'     \out{u<sub>ij</sub> = log2(k<sub>ij</sub>/s<sub>j</sub> + 1) and k<sub>ij</sub><sup>O</sup> = round(s<sub>j</sub>2<sup>μ<sub>i</sub><sup>u</sup>±exp(<i>&#119977</i>)σ<sub>i</sub><sup>u</sup></sup>)}}{
#'     \deqn{u_{ij}=log_2⁡(k_{ij}/s_{j} +1) and k_{ij}^{O}=round(s_{j} 2^{μ_{i}^{u}±exp(N) σ_{j}^{u}}}}.
#' More detail in the ABEILLE paper.
#' @importFrom matrixStats colSds
#' @examples
#' SequencingTable <- ExampleAbeilleDataSet
#' anno_table <- data.frame(sampleID = colnames(ExampleAbeilleDataSet),
#' Sex = sample(c("Male", "Female"), 128, replace = TRUE), Age = sample(1:89, 128, replace = TRUE),
#' Cohort = rep("GTEx database",128))
#' sizefactor <- SizeFactors(ExampleAbeilleDataSet, anno_table)
#' TableWithOutliers <- AddOutriderOutliers(SequencingTable, sizefactor, Zmean = log(3),
#' Zsd = log(1.6), drawoption = 10^-4, feedback = TRUE, setseed = 5)
#' print(head(TableWithOutliers$table))
#' print(head(TableWithOutliers$feedback))
#' @export
AddOutriderOutliers <- function(SequencingTable, SizeFactor, Zmean = log(3), Zsd = log(1.6), drawoption = 10^-5, feedback = FALSE, setseed = FALSE){
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
  uij <- log2(sweep(SequencingTable, 2, SizeFactor, FUN = '/') + 1)
  Znorm <- exp(rnorm(dim(SequencingTable)[1], Zmean, Zsd))
  stduij <- colSds(uij)
  randomsign <- matrix(sample(c(-1,1),dim(SequencingTable)[1]*dim(SequencingTable)[2],replace=TRUE),nrow=dim(SequencingTable)[1],ncol=dim(SequencingTable)[2])
  part1 <- randomsign * as.matrix(Znorm)%*%t(as.matrix(stduij))
  part2 <- 2**(uij + part1)
  onlyoutliers <- round(sweep(part2, 2, SizeFactor, FUN = '*'))
  original_value <- c()
  transformed_value <- c()
  for (i in 1:length(random_row)){
    original_value <- c(original_value,SequencingTable[random_row[i],random_col[i]])
    TableWithOutliers[random_row[i],random_col[i]] <- min(max(SequencingTable)+1,onlyoutliers[random_row[i],random_col[i]])
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
