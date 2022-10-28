#' Install missing packages, mainly for internal use.
#'
#' @return NULL
#' @export
#' @examples
#' require(mammalMethylationPredictors)
#' mammalMethylationPredictors::installMissingPackages()
#' @references
#' C. Li et al., Epigenetic predictors of maximum lifespan and
#' other life history traits in mammals. bioRxiv,  (2021).
#' A. Arneson et al., A mammalian methylation array for profiling methylation
#' levels at conserved sequences. Nature Communications 13, 1-13 (2022).



installMissingPackages <- function() {
  list.of.packages <- c("glmnet", "randomForest")
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages) > 0) install.packages(new.packages)
}
