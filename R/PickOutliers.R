#' Function to pick outliers
#'
#' @description Successively using linear regression and decision tree this function is able to pick outliers from two metric dataset. originally built to work a the z-score and the log2 fold-change it is possible to feed the function with others metrics. The decision tree can be also custom.
#'
#' @param zscore data frame of the same size as the initial dataset. The z-score can be compute with the function Zscore.
#' @param l2fc data frame of the same size as the initial dataset. The log2 fold-change can be compute with the function L2FC.
#' @param SequencingTable a read counts table with the transcripts in row and the samples in column.
#' @param ReconstructedTable the reconstructed table generate by an autoencoder.
#' @param decision_tree function containing a custom decision tree. Check the help of the function my_tree for more information.
#' @param nb_pval argument allowing to have an additional column in the result. The column is an estimation of the P-value using negative binomial distribution.
#'
#' @return Data frame containing only rows predict as outliers.
#' @importFrom stats cooks.distance
#' @importFrom stats dfbetas
#' @importFrom stats hatvalues
#' @importFrom stats lm
#' @importFrom stats rstandard
#' @importFrom stats complete.cases
#' @importFrom stats pnbinom
#' @importFrom progress progress_bar
#' @importFrom ABEILLE ComputeAberrantScore
#'
#' @examples
#' SequencingTable <- ExampleAbeilleDataSet
#' ReconstructedTable <- ExampleAbeilleReconstructed
#' zscore <- Zscore(SequencingTable,ReconstructedTable)
#' l2fc <- L2FC(SequencingTable,ReconstructedTable)
#' genedisp <- GeneDisp(SequencingTable)
#' outliers <- PickOutliers(zscore,l2fc,SequencingTable, ReconstructedTable,
#' decision_tree = my_tree, nb_pval=genedisp)
#' print(head(outliers))
#' @export




#' Function computing Aberrant Score based on the decision tree
#'
#' @description This function is call from the function PickOutliers
#'
#' @param zscore data frame of the same size as the initial dataset. The z-score can be compute with the function Zscore.
#'
#' @return Data frame with additional column containing Aberrant Score.
#'
#' @examples
#'
#' @export

ComputeAberrantScore <- function(my_df){

extra_col <- c()

for (i in seq(nrow(my_df))){

  if (my_df[i,"dfbetas_var"] >= 1.1){

    if(my_df[i,"hat"] >= 0.2){extra_col <- c(extra_col,4)}

    else{
      if (my_df[i,"typeerror"] < 2.4){
        extra_col <- c(extra_col, 1)
      }

      else{
        if (my_df[i,"l2fc"] >= 2.4){
          if (my_df[i,"zscore"] >= 6.6){
            extra_col <- c(extra_col, 2)
          }
        }
      }
    }

  }

  else{

    if(my_df[i,"zscore"] >= 11){
      extra_col <- c(extra_col, 3)
    }

    else if(my_df[i,"hat"] < 0.2 & my_df[i,"typeerror"] < 2.4){
      extra_col <- c(extra_col, 1)
    }

    else if(my_df[i,"hat"] < 0.2 & my_df[i,"typeerror"] >= 2.4 & my_df[i,"l2fc"] >= 2.4 & my_df[i,"zscore"] >= 6.6){
      extra_col <- c(extra_col, 2)
    }
  }

}
my_df$AberrantScore <- extra_col
return(my_df)
}

PickOutliers <- function(zscore, l2fc, SequencingTable, ReconstructedTable, decision_tree = my_tree, nb_pval = FALSE){
  if (nb_pval[1] == FALSE){
  i=1
  var1 <- zscore[i,]
  var2 <- l2fc[i,]
  linear_reg <- lm(var2~var1)
  DataReturn <- data.frame(Sample = NA,
                           Transcript = NA,
                           value = NA,
                           reconstruction = NA,
                           zscore=NA,
                           l2fc=NA,
                           typeerror=NA,
                           cooksD=NA,
                           hat=NA,
                           dfbetas_intercept=NA,
                           dfbetas_var=NA,
                           predict=NA)
  LinearRegression_data <- data.frame(Sample = colnames(SequencingTable),
                                      Transcript = rep(row.names(SequencingTable)[i],ncol(SequencingTable)),
                                      value = as.numeric(SequencingTable[i,]),
                                      reconstruction = as.numeric(ReconstructedTable[i,]),
                                      zscore=var1,
                                      l2fc=var2,
                                      typeerror=rstandard(linear_reg),
                                      cooksD=cooks.distance(linear_reg),
                                      hat=hatvalues(linear_reg),
                                      dfbetas_intercept=dfbetas(linear_reg)[,1],
                                      dfbetas_var=dfbetas(linear_reg)[,2],
                                      row.names = NULL)
  LinearRegression_data['predict'] <- decision_tree(LinearRegression_data)
  LinearRegression_data <- LinearRegression_data[LinearRegression_data$predict == TRUE,]
  LinearRegression_data <- LinearRegression_data[complete.cases(LinearRegression_data),]
  DataReturn <- rbind(DataReturn,LinearRegression_data)
  pb <- progress_bar$new(total = dim(SequencingTable)[1] - 1)
  pb$tick(0)
  for (i in 2:dim(SequencingTable)[1]){
    pb$tick()
    var1 <- zscore[i,]
    var2 <- l2fc[i,]
    linear_reg <- lm(var2~var1)
    LinearRegression_data <- data.frame(Sample = colnames(SequencingTable),
                                        Transcript = rep(row.names(SequencingTable)[i],ncol(SequencingTable)),
                                        value = as.numeric(SequencingTable[i,]),
                                        reconstruction = as.numeric(ReconstructedTable[i,]),
                                        zscore=var1,
                                        l2fc=var2,
                                        typeerror=rstandard(linear_reg),
                                        cooksD=cooks.distance(linear_reg),
                                        hat=hatvalues(linear_reg),
                                        dfbetas_intercept=dfbetas(linear_reg)[,1],
                                        dfbetas_var=dfbetas(linear_reg)[,2],
                                        row.names = NULL)
    LinearRegression_data['predict'] <- decision_tree(LinearRegression_data)
    LinearRegression_data <- LinearRegression_data[LinearRegression_data$predict == TRUE,]
    LinearRegression_data <- LinearRegression_data[complete.cases(LinearRegression_data), ]
    DataReturn <- rbind(DataReturn,LinearRegression_data)}
  DataReturn <- DataReturn[-1,]
  row.names(DataReturn) <- NULL
  DataReturn <- ComputeAberrantScore(DataReturn)
  return(DataReturn)
  }else if(class(nb_pval[1]) == "numeric"){
    if (class(SequencingTable) == "data.frame"){SequencingTable <- as.matrix(SequencingTable)}
    if (class(ReconstructedTable) == "data.frame"){ReconstructedTable <- as.matrix(ReconstructedTable)}
    i=1
    var1 <- zscore[i,]
    var2 <- l2fc[i,]
    linear_reg <- lm(var2~var1)
    DataReturn <- data.frame(Sample = NA,
                             Transcript = NA,
                             value = NA,
                             reconstruction = NA,
                             zscore=NA,
                             l2fc=NA,
                             typeerror=NA,
                             cooksD=NA,
                             hat=NA,
                             dfbetas_intercept=NA,
                             dfbetas_var=NA,
                             nb_pval = NA,
                             predict=NA)
    LinearRegression_data <- data.frame(Sample = colnames(SequencingTable),
                                        Transcript = rep(row.names(SequencingTable)[i],ncol(SequencingTable)),
                                        value = as.numeric(SequencingTable[i,]),
                                        reconstruction = as.numeric(ReconstructedTable[i,]),
                                        zscore=var1,
                                        l2fc=var2,
                                        typeerror=rstandard(linear_reg),
                                        cooksD=cooks.distance(linear_reg),
                                        hat=hatvalues(linear_reg),
                                        dfbetas_intercept=dfbetas(linear_reg)[,1],
                                        dfbetas_var=dfbetas(linear_reg)[,2],
                                        nb_pval = pnbinom(q = SequencingTable[i,], size = nb_pval[i], mu = ReconstructedTable[i,]),
                                        row.names = NULL)
    LinearRegression_data['predict'] <- decision_tree(LinearRegression_data)
    LinearRegression_data <- LinearRegression_data[LinearRegression_data$predict == TRUE,]
    LinearRegression_data <- LinearRegression_data[complete.cases(LinearRegression_data),]
    DataReturn <- rbind(DataReturn,LinearRegression_data)
    pb <- progress_bar$new(total = dim(SequencingTable)[1] - 1)
    pb$tick(0)
    for (i in 2:dim(SequencingTable)[1]){
      pb$tick()
      var1 <- zscore[i,]
      var2 <- l2fc[i,]
      linear_reg <- lm(var2~var1)
      LinearRegression_data <- data.frame(Sample = colnames(SequencingTable),
                                          Transcript = rep(row.names(SequencingTable)[i],ncol(SequencingTable)),
                                          value = as.numeric(SequencingTable[i,]),
                                          reconstruction = as.numeric(ReconstructedTable[i,]),
                                          zscore=var1,
                                          l2fc=var2,
                                          typeerror=rstandard(linear_reg),
                                          cooksD=cooks.distance(linear_reg),
                                          hat=hatvalues(linear_reg),
                                          dfbetas_intercept=dfbetas(linear_reg)[,1],
                                          dfbetas_var=dfbetas(linear_reg)[,2],
                                          nb_pval = pnbinom(q = SequencingTable[i,], size = nb_pval[i], mu = ReconstructedTable[i,]),
                                          row.names = NULL)
      LinearRegression_data['predict'] <- decision_tree(LinearRegression_data)
      LinearRegression_data <- LinearRegression_data[LinearRegression_data$predict == TRUE,]
      LinearRegression_data <- LinearRegression_data[complete.cases(LinearRegression_data), ]
      DataReturn <- rbind(DataReturn,LinearRegression_data)}
    DataReturn <- DataReturn[-1,]
    row.names(DataReturn) <- NULL
    DataReturn <- ComputeAberrantScore(DataReturn)
    return(DataReturn)
  }else{
    stop("Numerical vector or logical needed in nb_pval")
  }
}
