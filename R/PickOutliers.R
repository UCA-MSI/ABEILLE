#' @title Function to pick outliers
#'
#' @description Successively using linear regression and decision tree this function is able to pick outliers from two metric dataset. originally built to work a the z-score and the log2 fold-change it is possible to feed the function with others metrics. The decision tree can be also custom.
#'
#' @param divergence_score data frame of the same size as the initial dataset. The divergence score can be compute with the function divergence_score.
#' @param delta_count data frame of the same size as the initial dataset. The log2 fold-change can be compute with the function delta_count.
#' @param SequencingTable a read counts table with the transcripts in row and the samples in column.
#' @param ReconstructedTable the reconstructed table generate by an autoencoder.
#' @param decision_tree function containing a custom decision tree. Check the help of the function my_tree for more information.
#'
#' @return Data frame containing only rows predict as outliers.
#' @importFrom stats cooks.distance
#' @importFrom stats dfbetas
#' @importFrom stats hatvalues
#' @importFrom stats lm
#' @importFrom stats rstandard
#' @importFrom stats complete.cases
#' @importFrom progress progress_bar
#'
#' @examples
#' SequencingTable <- ExampleAbeilleDataSet
#' ReconstructedTable <- ExampleAbeilleReconstructed
#' divergence_score <- DivergenceScore(SequencingTable,ReconstructedTable)
#' delta_count <- DeltaCount(SequencingTable,ReconstructedTable)
#' genedisp <- GeneDisp(SequencingTable)
#' outliers <- PickOutliers(divergence_score,delta_count,SequencingTable, ReconstructedTable,
#' decision_tree = my_tree)
#' print(head(outliers))
#' @export

PickOutliers <- function(divergence_score, delta_count, SequencingTable, ReconstructedTable, decision_tree = my_tree){
  i=1
  var1 <- divergence_score[i,]
  var2 <- delta_count[i,]
  linear_reg <- lm(var2~var1)
  DataReturn <- data.frame(Sample = NA,
                           Transcript = NA,
                           value = NA,
                           reconstruction = NA,
                           divergence_score=NA,
                           delta_count=NA,
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
                                      divergence_score=var1,
                                      delta_count=var2,
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
    var1 <- divergence_score[i,]
    var2 <- delta_count[i,]
    linear_reg <- lm(var2~var1)
    LinearRegression_data <- data.frame(Sample = colnames(SequencingTable),
                                        Transcript = rep(row.names(SequencingTable)[i],ncol(SequencingTable)),
                                        value = as.numeric(SequencingTable[i,]),
                                        reconstruction = as.numeric(ReconstructedTable[i,]),
                                        divergence_score=var1,
                                        delta_count=var2,
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
}

#' @section
#'
#' @title Function computing Aberrant Score based on the decision tree
#'
#' @name ComputeAberrantScore
#' @description This function is call from the function PickOutliers
#'
#' @param divergence_score data frame of the same size as the initial dataset. The z-score can be compute with the function divergence_score.
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
        if (my_df[i,"delta_count"] >= 2.4){
          if (my_df[i,"divergence_score"] >= 6.6){
            extra_col <- c(extra_col, 2)
          }
        }
      }
    }

  }

  else{

    if(my_df[i,"divergence_score"] >= 11){
      extra_col <- c(extra_col, 3)
    }

    else if(my_df[i,"hat"] < 0.2 & my_df[i,"typeerror"] < 2.4){
      extra_col <- c(extra_col, 1)
    }

    else if(my_df[i,"hat"] < 0.2 & my_df[i,"typeerror"] >= 2.4 & my_df[i,"delta_count"] >= 2.4 & my_df[i,"divergence_score"] >= 6.6){
      extra_col <- c(extra_col, 2)
    }
  }

}
my_df$AberrantScore <- extra_col
return(my_df)
}