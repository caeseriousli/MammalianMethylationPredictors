#' Predict lifespan, given normalized Mammalian array data.
#'
#' @param dt an nxp matrix, or data frame; normalized Mammalian array data; must have been normalized by SeSAme normalization pipeline
#' @param version specifies which variant of the predictor to us; default to "40K"; can be one of c("40K", "320K40K")
#' @param predictorType one of c("final", "EPIC") return just a vector or a matrix including second guess and their prediction probabilities
#' @return a vector or a matrix of predicted values
#' @export
#' @examples
#' require(mammalMethylationPredictors)
#' # Load your normalized DNA methylation data. For normalization pipelines, see shorvath/MammalianMethylationConsortium (Arneson)
#' dat0 <- readRDS("PATH_TO_YOUR_DATA")
#'
#' # fit the predictor
#' results = predictLifespan(dt = dat0, arrayType = "40K")
#'
#' @references
#' A. Arneson et al., A mammalian methylation array for profiling methylation
#' levels at conserved sequences. Nature Communications 13, 1-13 (2022).
#' C. Li et al., Epigenetic predictors of maximum lifespan and
#' other life history traits in mammals. bioRxiv,  (2021).

predictLifespan <- function(dt = NULL, arrayType = "40K", predictorType = "final") {

    if(arrayType == "40K" & ncol(dt) > 4e4) {
      arrayType = "320K"
    } else if(arrayType == "320K" & ncol(dt) < 1e5) {
      arrayType = "40K"
    }



    # Progress bar
    pb = txtProgressBar(min = 0, max = 3, initial = 0, style = 3)
    setTxtProgressBar(pb,1)

    # sort(sapply(fit$fit,function(x){object.size(x)}))

    ## Load the Sex Predictor
    fit = read.csv(paste0("./Predictors/LifespanPredictor_40K_Li2021.csv"), stringsAsFactors = FALSE)

    # Progress bar
    setTxtProgressBar(pb,2)

    if(arrayType == "40K") {
      # dt = dt[, fit[-1, "CpG"]]

    } else if(arrayType == "320K") {
      dictionary = readRDS("./inst/Mapping_OneToOneProbes_320K_40K.RDS")

      ## Note that the feature names used in 1-1 fit$featureNames ensures CGid are unique
      ## Now First re-order the Amin dictionary to translate 320K colnames to RF feature names (40K CGid)
      # if(sum(!colnames(dt) %in% dictionary$Probe_ID) > 0) stop("Some probes not in the dictionary")
      rownames(dictionary) = dictionary$Probe_ID; dictionary = dictionary[colnames(dt), ]
      colnames(dt) = dictionary$CGid

      ## Imputing dat0 for a few 40K probes that do not exist on 320K
      # Reduce dt's length first
      dt = dt[, colnames(dt) %in% fit[, "CpG"]]
      #
      myimpute = fit[-1, "CpG"]; myimpute = myimpute[!myimpute %in% colnames(dt)]
      myimpute_matrix = matrix(0.5, nrow = nrow(dt), ncol = length(myimpute)); colnames(myimpute_matrix) = myimpute
      dt = cbind(dt, myimpute_matrix)
      ######

    }

    if(sum(!fit[-1, "CpG"] %in% colnames(dt)) > 0) {
      stop("Double check arrayType argument or data. The arrayType probes do not match our predictor.")
    }

    dt = dt[, fit[-1, "CpG"]]

    # Progress bar
    setTxtProgressBar(pb,3)
    close(pb)

    return(apply(dt, 1, function(x) return(sum(x * fit[-1, 1]) + fit[1, 1])))

}

