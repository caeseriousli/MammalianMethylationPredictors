#' Predict Tissue, given normalized Mammalian array data.
#'
#' @param dt an nxp matrix, or data frame; normalized Mammalian array data; must have been normalized by SeSAme normalization pipeline
#' @param arrayType specifies which variant of the predictor to us; default to "40K"; can be one of c("40K", "320K")
#' @param returnType one of c("vector", "moreInfo") return just a vector or a matrix including second guess and their prediction probabilities
#' @return a vector or a matrix of predicted values
#' @export
#' @examples
#' require(mammalMethylationPredictors)
#' # Load your normalized DNA methylation data. For normalization pipelines, see shorvath/MammalianMethylationConsortium (Arneson)
#' mydata <- readRDS("PATH_TO_YOUR_DATA")
#'
#' # fit the predictor
#' results = mammalMethylationPredictors::predictTissue(dt = mydata, arrayType = "40K", returnType = "moreInfo")
#'
#' @references
#' C. Li et al., Epigenetic predictors of maximum lifespan and other life history traits in mammals. bioRxiv, (2021).
#'
#' Haghani, et al., DNA Methylation Networks Underlying Mammalian Traits, bioRxiv, (2021).
#'
#' A. T. Lu et al., Universal DNA methylation age across mammalian tissues. bioRxiv, 2021.2001.2018.426733 (2021)
#'
#' A. Arneson et al., A mammalian methylation array for profiling methylation levels at conserved sequences. Nature Communications 13, 1-13 (2022).

predictTissue <- function(dt = NULL, arrayType = "40K", returnType = "moreInfo") {

  installMissingPackages()
  library(randomForest)
  # Progress bar
  pb = txtProgressBar(min = 0, max = 3, initial = 0, style = 3)
  setTxtProgressBar(pb,1)

  predictWhat = "Tissue"
  if(arrayType == "40K" & ncol(dt) > 4e4) {
    stop("Double check arrayType argument. It is specified as 40K array but dt has more than 40 thousand probes.")
  }

  # sort(sapply(fit$fit,function(x){object.size(x)}))

  # Progress
  setTxtProgressBar(pb,2)

  ## Load the Tissue Predictor
  mydata <- system.file("extdata", "Tissue_Overlap320K40K_Filter5samples_100trees_randomForest.RDS", package = "mammalMethylationPredictors")
  fit = readRDS(mydata)
  if(arrayType == "40K") {
    dt = dt[, fit$featureNames]

  } else if(arrayType == "320K") {
    mydata <- system.file("extdata", "Mapping_OneToOneProbes_320K_40K.RDS", package = "mammalMethylationPredictors")
    dictionary = readRDS(mydata)

    ## Note that the feature names used in 1-1 fit$featureNames ensures CGid are unique
    ## Now First re-order the Amin dictionary to translate 320K colnames to RF feature names (40K CGid)
    # if(sum(!colnames(dt) %in% dictionary$Probe_ID) > 0) stop("Some probes not in the dictionary")
    rownames(dictionary) = dictionary$Probe_ID; dictionary = dictionary[colnames(dt), ]
    colnames(dt) = dictionary$CGid

    dt = dt[, fit$featureNames]
  }

  # Progress
  setTxtProgressBar(pb,3)
  close(pb)

  if(returnType == "moreInfo") {
    temp = predict(fit$fit, dt, type="prob")
    samples = data.frame(matrix(NA, ncol = 4, nrow = nrow(dt)))
    colnames(samples) = c(paste0("Predicted_", predictWhat), paste0("Predicted_", predictWhat, "_Probability"),
                          paste0("Predicted_", predictWhat, "_SecondGuess"), paste0("Predicted_", predictWhat, "_SecondProbability"))

    samples[, paste0("Predicted_", predictWhat)] <- predict(fit$fit, dt, type="response")

    samples[, paste0("Predicted_", predictWhat, "_Probability")] <-
      apply(temp, 1, function(x) {
        return(max(x))
      })

    samples[, paste0("Predicted_", predictWhat, "_SecondGuess")] = apply(temp, 1, function(x) {
      allclasses = colnames(temp)
      allclasses = allclasses[order(x, decreasing = T)]
      return(allclasses[2])
    })

    samples[, paste0("Predicted_", predictWhat, "_SecondProbability")] = apply(temp, 1, function(x) {
      x = x[order(x, decreasing = T)]
      return(x[2])
    })

    return(samples)
  } else {
    return(predict(fit$fit, dt, ty))
  }

}

