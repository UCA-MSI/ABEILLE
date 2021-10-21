#' Decision tree to choose outliers
#'
#' @description The decision tree to pick outliers is a function using filters on different variables of the linear regression. This one was created by using CART on OUTRIDER outliers on the GTEx database. It is possible to use a custom decision tree by chaining ifelse conditions. Here are all variables usable:
#' \itemize{
#'  \item{"value"}{ the original value.}
#'  \item{"reconstruction"}{ the value reconstructed.}
#'  \item{"divergence_score"}{ see DivergenceScore function.}
#'  \item{"delta_count"}{ see DeltaCount function.}
#'  \item{"typeerror"}{ standardized residual. For more check rstandard function.}
#'  \item{"cooksD"}{ cook's distance. For more check cooks.distance function.}
#'  \item{"hat"}{ projection matrix. For more check hatvalues function.}
#'  \item{"dfbetas_intercept and dfbetas_var"}{ influence measures on the linear regression. For more check dfbetas function.}
#' }
#'
#' @param linear_reg_data data of the linear regression.
#'
#' @return Boolean vector. This function is build to feed the PickOutliers one, but you can build your own.
#' @examples
#' \dontrun{
#'  #The decision tree used:
#'  bool_vector <- ifelse(linear_reg_data['dfbetas_var'] >= 1.1,
#'  ifelse(linear_reg_data['hat'] >= 0.2, TRUE,
#'  ifelse(linear_reg_data['typeerror'] < 2.4, TRUE,
#'  ifelse(linear_reg_data['delta_count'] >= 2.4,
#'  ifelse(linear_reg_data['divergence_score'] >= 6.6, TRUE, FALSE), FALSE))),
#'  ifelse(linear_reg_data['divergence_score'] >=11, TRUE, FALSE))
#'
#'  #A simple tree possible:
#'  bool_vector <- ifelse(linear_reg_data['cooksD'] >= 5, TRUE, FALSE)
#' }
#' @export

my_tree <- function(linear_reg_data){
  bool_vector <- ifelse(linear_reg_data['dfbetas_var'] >= 1.1,
                        ifelse(linear_reg_data['hat'] >= 0.2, TRUE, ifelse(linear_reg_data['typeerror'] < 2.4, TRUE, ifelse(linear_reg_data['delta_count'] >= 2.4, ifelse(linear_reg_data['divergence_score'] >= 6.6, TRUE, FALSE), FALSE))),
                        ifelse(linear_reg_data['divergence_score'] >=11,TRUE,FALSE))
  return(bool_vector)
}
