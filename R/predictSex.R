#' Predict Sex, given normalized Mammalian array data.
#'
#' @param dt an nxp matrix, or data frame; normalized Mammalian array data; must have been normalized by SeSAme normalization pipeline
#' @param arrayType specifies which variant of the predictor to us; default to "40K"; can be one of c("40K", "320K")
#' @param returnType one of c("vector", "moreInfo") return just a vector or a matrix including second guess and their prediction probabilities
#' @return a vector or a matrix of predicted values
#' @export
#' @examples
#' require(mammalMethylationPredictors)
#' # Load your normalized DNA methylation data. For normalization pipelines, see shorvath/MammalianMethylationConsortium (Arneson, 2021)
#' mydata <- readRDS("PATH_TO_YOUR_DATA")
#'
#' # fit the predictor
#' results = predictSex(dt = mydata, arrayType = "40K", returnType = "moreInfo")
#'
#' @references
#' C. Li et al., Epigenetic predictors of maximum lifespan and other life history traits in mammals. bioRxiv, (2021).
#'
#' Haghani, et al., DNA Methylation Networks Underlying Mammalian Traits, bioRxiv, (2021).
#'
#' A. T. Lu et al., Universal DNA methylation age across mammalian tissues. bioRxiv, 2021.2001.2018.426733 (2021)
#'
#' A. Arneson et al., A mammalian methylation array for profiling methylation levels at conserved sequences. Nature Communications 13, 1-13 (2022).


predictSex <- function(dt = NULL, arrayType = "40K", returnType = "moreInfo") {

  if(arrayType == "40K" & ncol(dt) > 4e4) {
    stop("Double check arrayType argument. It is specified as 40K array but dt has more than 40 thousand probes.")
  }

  # Progress bar
  pb = txtProgressBar(min = 0, max = 3, initial = 0, style = 3)
  setTxtProgressBar(pb,1)

  # sort(sapply(fit$fit,function(x){object.size(x)}))

  ## Load the Sex Predictor
  mydata <- system.file("extdata", "FemalePredictor_Overlap320K40K.csv", package = "mammalMethylationPredictors")
  fit = read.csv(mydata, stringsAsFactors = FALSE)
  # fit = read.csv(paste0("./Predictors/FemalePredictor_Overlap320K40K.csv"), stringsAsFactors = FALSE)

  # Progress bar
  setTxtProgressBar(pb,2)

  if(arrayType == "40K") {
    dt = dt[, fit[-1, "CpG"]]

  } else if(arrayType == "320K") {
    mydata <- system.file("extdata", "Mapping_OneToOneProbes_320K_40K.RDS", package = "mammalMethylationPredictors")
    dictionary = readRDS(mydata)

    ## Note that the feature names used in 1-1 fit$featureNames ensures CGid are unique
    ## Now First re-order the Amin dictionary to translate 320K colnames to RF feature names (40K CGid)
    # if(sum(!colnames(dt) %in% dictionary$Probe_ID) > 0) stop("Some probes not in the dictionary")
    rownames(dictionary) = dictionary$Probe_ID; dictionary = dictionary[colnames(dt), ]
    colnames(dt) = dictionary$CGid

    dt = dt[, fit[-1, "CpG"]]
  }

  # Progress bar
  setTxtProgressBar(pb,3)
  close(pb)

  if(returnType == "moreInfo") {

    samples = data.frame(matrix(NA, ncol = 2, nrow = nrow(dt)))
    colnames(samples) = c(paste0("Predicted_", "Female"),
                          paste0("Predicted_", "Female_Probability"))

    temp = apply(dt, 1, function(x) return(sum(x * fit[-1, 1]) + fit[1, 1]))
    temp = exp(temp) / (1+exp(temp))
    samples[, paste0("Predicted_", "Female")] = ifelse(temp >= 0.5, 1, 0)
    samples[, paste0("Predicted_", "Female_Probability")] <- temp

    return(samples)
  } else {
    temp = apply(dt, 1, function(x) return(sum(x * fit[-1, 1]) + fit[1, 1]))
    temp = exp(temp) / (1+exp(temp))
    return(ifelse(temp >= 0.5, 1, 0))
  }



}

