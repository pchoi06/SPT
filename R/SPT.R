#' @title Sigmoidal vs. Double-Sigmoidal Model Selection Tool
#'
#' @description
#' Fit both a single-sigmoidal and a double-sigmoidal curve to time-intensity data,
#' then choose which model best describes the data—or report no signal or ambiguous.
#'
#' @param dataInput A data frame (columns `time`, `intensity`) or a list with
#'   element `timeIntensityData`.
#' @param dataInputName Optional character name for the dataset, used in the output.
#' @param n_runs_min_sm Minimum number of successful sigmoidal fits before stopping.
#' @param n_runs_max_sm Maximum total sigmoidal fit attempts.
#' @param n_runs_min_dsm Minimum number of successful double-sigmoidal fits before stopping.
#' @param n_runs_max_dsm Maximum total double-sigmoidal fit attempts.
#' @param showDetails Logical; if TRUE, print intermediate objects for debugging.
#' @param startList_sm,lowerBounds_sm,upperBounds_sm,min_Factor_sm,n_iterations_sm
#'   Sigmoidal fit control: initial guesses, optimizer bounds, minimal step factor,
#'   and maximum iterations (passed to `sigFit`).
#' @param startList_dsm,lowerBounds_dsm,upperBounds_dsm,min_Factor_dsm,n_iterations_dsm
#'   Double-sigmoidal fit control: analogous parameters for `sigFit2`.
#' @param threshold_intensity_range Numeric; minimum data dynamic range required.
#' @param threshold_minimum_for_intensity_maximum Numeric; minimum peak height required.
#' @param threshold_bonus_sigmoidal_AIC Numeric; AIC offset favoring sigmoidal.
#' @param threshold_sm_tmax_IntensityRatio Numeric; ratio of sigmoidal intensity at max time.
#' @param threshold_dsm_tmax_IntensityRatio Numeric; ratio of double-sigmoidal intensity at max time.
#' @param threshold_AIC Numeric; minimum AIC improvement to accept a fit.
#' @param threshold_t0_max_int Numeric; maximum allowed baseline intensity.
#' @param stepSize Numeric; finite‐difference step for slope calculations.
#' @param ... Additional arguments passed to underlying fit functions.
#'
#' @return A list with elements:
#' \describe{
#'   \item{normalizedInput}{The normalized input data (list).}
#'   \item{sigmoidalModel}{Data frame of sigmoidal fit parameters.}
#'   \item{doubleSigmoidalModel}{Data frame of double-sigmoidal fit parameters.}
#'   \item{decisionProcess}{List of logical tests and thresholds used.}
#'   \item{summaryVector}{Named vector summarizing the chosen model and key parameters.}
#' }
#' @export
SPT <- function (dataInput, dataInputName = NA, n_runs_min_sm = 20,
                 n_runs_max_sm = 500, n_runs_min_dsm = 20, n_runs_max_dsm = 500,
                 showDetails = FALSE,
                 startList_sm = list(h0 = 0, maximum = 1, slopeParam = 1, midPoint = 0.33),
                 lowerBounds_sm = c(h0 = 0, maximum = 0.3, slopeParam = 0.01, midPoint = -0.52),
                 upperBounds_sm = c(h0 = 0.3, maximum = 1.5, slopeParam = 180,  midPoint = 1.15),
                 min_Factor_sm = 1/2^20, n_iterations_sm = 1000,
                 startList_dsm = list(
                   finalAsymptoteIntensityRatio = 0.5,
                   maximum                      = 1,
                   slope1Param                  = 1,
                   midPoint1Param               = 0.33,
                   slope2Param                  = 1,
                   midPointDistanceParam        = 0.29,
                   h0                           = 0.1     # ← added
                 ),
                 lowerBounds_dsm = c(
                   finalAsymptoteIntensityRatio = 0,
                   maximum                      = 0.3,
                   slope1Param                  = 0.01,
                   midPoint1Param               = -0.52,
                   slope2Param                  = 0.01,
                   midPointDistanceParam        = 0.04,
                   h0                           = 0    # ← added
                 ),
                 upperBounds_dsm = c(
                   finalAsymptoteIntensityRatio = 1.5,
                   maximum                      = 1.5,
                   slope1Param                  = 180,
                   midPoint1Param               = 1.15,
                   slope2Param                  = 180,
                   midPointDistanceParam        = 0.63,
                   h0                           = 0.3    # ← added
                 ),
                 min_Factor_dsm = 1/2^20, n_iterations_dsm = 1000, threshold_intensity_range = 0.1,
                 threshold_minimum_for_intensity_maximum = 0.3, threshold_bonus_sigmoidal_AIC = 0,
                 threshold_sm_tmax_IntensityRatio = 0.85, threshold_dsm_tmax_IntensityRatio = 0.75,
                 threshold_AIC = -10, threshold_t0_max_int = 0.05, stepSize = 1e-05,
                 ...)
{
  normalizedInput = normalizeData_spt(dataInput = dataInput,
                                           dataInputName = dataInputName)
  preDecisionProcess = preCategorize_spt(normalizedInput = normalizedInput,
                                              threshold_intensity_range = threshold_intensity_range,
                                              threshold_minimum_for_intensity_maximum = threshold_minimum_for_intensity_maximum)
  if (showDetails) {
    utils::str(preDecisionProcess)
  }
  if (preDecisionProcess$decision == "no_signal") {
    summaryVector <- c()
    if (!is.na(normalizedInput$dataInputName)) {
      summaryVector$dataInputName <- as.vector(normalizedInput$dataInputName)
    }
    else {
      summaryVector$dataInputName <- NA
    }
    summaryVector$decision <- "no_signal"
    return(list(normalizedInput = normalizedInput, sigmoidalModel = NA,
                doubleSigmoidalModel = NA, decisionProcess = preDecisionProcess,
                summaryVector = summaryVector))
  }
  if (preDecisionProcess$decision == "not_no_signal") {
    sigmoidalModel <- multipleFitFunction_spt(dataInput = normalizedInput,
                                  model = "sigmoidal", n_runs_min = n_runs_min_sm,
                                  n_runs_max = n_runs_max_sm, showDetails = showDetails,
                                  startList = startList_sm, lowerBounds = lowerBounds_sm,
                                  upperBounds = upperBounds_sm, min_Factor = min_Factor_sm,
                                  n_iterations = n_iterations_sm)
    doubleSigmoidalModel <- multipleFitFunction_spt(dataInput = normalizedInput,
                                        model = "doublesigmoidal", n_runs_min = n_runs_min_dsm,
                                        n_runs_max = n_runs_max_dsm, showDetails = showDetails,
                                        startList = startList_dsm, lowerBounds = lowerBounds_dsm,
                                        upperBounds = upperBounds_dsm, min_Factor = min_Factor_dsm,
                                        n_iterations = n_iterations_dsm)
    sigmoidalModel <- parameterCalculation_spt(parameterVector = sigmoidalModel,
                                              stepSize = stepSize)
    doubleSigmoidalModel <- parameterCalculation_spt(parameterVector = doubleSigmoidalModel,
                                                    stepSize = stepSize)
    decisionProcess <- Categorize_spt(parameterVectorSigmoidal = sigmoidalModel,
                                     parameterVectorDoubleSigmoidal = doubleSigmoidalModel,
                                     threshold_intensity_range = threshold_intensity_range,
                                     threshold_minimum_for_intensity_maximum = threshold_minimum_for_intensity_maximum,
                                     threshold_bonus_sigmoidal_AIC = threshold_bonus_sigmoidal_AIC,
                                     threshold_sm_tmax_IntensityRatio = threshold_sm_tmax_IntensityRatio,
                                     threshold_dsm_tmax_IntensityRatio = threshold_dsm_tmax_IntensityRatio,
                                     threshold_AIC = threshold_AIC, threshold_t0_max_int = threshold_t0_max_int,
                                     showDetails = showDetails)
    summaryVector <- c()
    if (decisionProcess$decision == "sigmoidal") {
      summaryVector$dataInputName <- decisionProcess$dataInputName
      summaryVector$decision <- "sigmoidal"
      summaryVector$maximum_x <- sigmoidalModel$maximum_x
      summaryVector$maximum_y <- sigmoidalModel$maximum_y
      summaryVector$midPoint_x <- sigmoidalModel$midPoint_x
      summaryVector$midPoint_y <- sigmoidalModel$midPoint_y
      summaryVector$slope <- sigmoidalModel$slope
      summaryVector$incrementTime <- sigmoidalModel$incrementTime
      summaryVector$startPoint_x <- sigmoidalModel$startPoint_x
      summaryVector$startPoint_y <- sigmoidalModel$startPoint_y
      summaryVector$reachMaximum_x <- sigmoidalModel$reachMaximum_x
      summaryVector$reachMaximum_y <- sigmoidalModel$reachMaximum_y
    }
    if (decisionProcess$decision == "double_sigmoidal") {
      summaryVector$dataInputName <- decisionProcess$dataInputName
      summaryVector$decision <- "double_sigmoidal"
      summaryVector$maximum_x <- doubleSigmoidalModel$maximum_x
      summaryVector$maximum_y <- doubleSigmoidalModel$maximum_y
      summaryVector$midPoint1_x <- doubleSigmoidalModel$midPoint1_x
      summaryVector$midPoint1_y <- doubleSigmoidalModel$midPoint1_y
      summaryVector$midPoint2_x <- doubleSigmoidalModel$midPoint2_x
      summaryVector$midPoint2_y <- doubleSigmoidalModel$midPoint2_y
      summaryVector$slope1 <- doubleSigmoidalModel$slope1
      summaryVector$slope2 <- doubleSigmoidalModel$slope2
      summaryVector$finalAsymptoteIntensity <- doubleSigmoidalModel$finalAsymptoteIntensity
      summaryVector$incrementTime <- doubleSigmoidalModel$incrementTime
      summaryVector$startPoint_x <- doubleSigmoidalModel$startPoint_x
      summaryVector$startPoint_y <- doubleSigmoidalModel$startPoint_y
      summaryVector$reachMaximum_x <- doubleSigmoidalModel$reachMaximum_x
      summaryVector$reachMaximum_y <- doubleSigmoidalModel$reachMaximum_y
      summaryVector$decrementTime <- doubleSigmoidalModel$decrementTime
      summaryVector$startDeclinePoint_x <- doubleSigmoidalModel$startDeclinePoint_x
      summaryVector$startDeclinePoint_y <- doubleSigmoidalModel$startDeclinePoint_y
      summaryVector$endDeclinePoint_x <- doubleSigmoidalModel$endDeclinePoint_x
      summaryVector$endDeclinePoint_y <- doubleSigmoidalModel$endDeclinePoint_y
    }
    if (decisionProcess$decision == "no_signal") {
      summaryVector$dataInputName <- decisionProcess$dataInputName
      summaryVector$decision <- "no_signal"
    }
    if (decisionProcess$decision == "ambiguous") {
      summaryVector$dataInputName <- decisionProcess$dataInputName
      summaryVector$decision <- "ambiguous"
    }
    if (showDetails) {
      utils::str(decisionProcess)
    }
    return(list(normalizedInput = normalizedInput, sigmoidalModel = sigmoidalModel,
                doubleSigmoidalModel = doubleSigmoidalModel, decisionProcess = decisionProcess,
                summaryVector = summaryVector))
  }
}
