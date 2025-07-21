#' Plot Observed Data and Fitted Sigmoidal / Double-Sigmoidal Curves
#'
#' @description
#' Given raw time–intensity data and one or both fitted model parameter vectors,
#' plots the data points together with the theoretical sigmoidal and/or
#' double-sigmoidal curves. Optionally adds dashed lines or markers to
#' illustrate key model parameters (e.g. midpoints, plateaus).
#'
#' @param dataInput
#'   A data frame with columns `time` and `intensity`, or a list
#'   returned by the SPT fitting workflow (containing
#'   `$timeIntensityData` and scaling parameters).
#' @param sigmoidalFitVector
#'   *(Optional)* A fitted sigmoidal model vector produced by SPT (must have
#'   `model == "sigmoidal"`).
#' @param doubleSigmoidalFitVector
#'   *(Optional)* A fitted double-sigmoidal model vector produced by SPT
#'   (must have `model == "doublesigmoidal"`).
#' @param showParameterRelatedLines
#'   Logical. If `TRUE`, overlays horizontal/vertical segments and points
#'   that illustrate estimated parameters (requires that
#'   `additionalParameters == TRUE` in the fit vectors).
#' @param xlabelText
#'   Character. Label for the x-axis (default `"time"`).
#' @param ylabelText
#'   Character. Label for the y-axis (default `"intensity"`).
#' @param fittedXmin
#'   Numeric. Minimum x-value to use when drawing the theoretical curves;
#'   defaults to `0`.
#' @param fittedXmax
#'   Numeric or `NA`. Maximum x-value for the theoretical curves; if `NA`,
#'   uses the full time range from the fit vector.
#'
#' @return
#' A `ggplot` object showing the points and fitted curve(s).
#'
#' @examples
#' \dontrun{
#'   # simulate and fit, then plot:
#'   df <- tibble(time = seq(0, 300, 50),
#'                intensity = sigmoidalFitFormula_h0(seq(0,300,50),
#'                      maximum=480, slopeParam=0.04, midPoint=160, h0=0)
#'                         + rnorm(7,0,50))
#'   fit <- SPT(dataInput = df)
#'   figureModelCurves_h0(df, sigmoidalFitVector = fit$sigmoidalModel,
#'                         showParameterRelatedLines = TRUE)
#' }
#'
#' @import ggplot2
#' @export
figureModelCurves_h0  <- function (dataInput, sigmoidalFitVector = NULL, doubleSigmoidalFitVector = NULL,
                                  showParameterRelatedLines = FALSE, xlabelText = "time",
                                  ylabelText = "intensity", fittedXmin = 0, fittedXmax = NA)
{
  dataOutputVariable <- dataCheck_h0(dataInput)
  sameSourceDataCheck_h0(dataInput, sigmoidalFitVector, doubleSigmoidalFitVector)
  isalist <- (is.list(dataInput) & !is.data.frame(dataInput))
  if (isalist) {
    dataInput <- unnormalizeData_h0(dataInput)
    dataFrameInput <- dataInput$timeIntensityData
  }
  isadataframe = (is.data.frame(dataInput))
  if (isadataframe) {
    dataFrameInput <- dataInput
  }
  if (!is.null(sigmoidalFitVector)) {
    if (!sigmoidalFitVector$model == "sigmoidal") {
      stop("provided sigmoidalFitVector is not a sigmoidal fit vector")
    }
    if (!sigmoidalFitVector$isThisaFit) {
      warning("provided sigmoidal fit vector does not include a fit!")
    }
    if (sigmoidalFitVector$isThisaFit) {
      maximum_x <- sigmoidalFitVector$maximum_x
      maximum_y <- sigmoidalFitVector$maximum_y
      midPoint_x <- sigmoidalFitVector$midPoint_x
      midPoint_y <- sigmoidalFitVector$midPoint_y
      slope <- sigmoidalFitVector$slope
      slopeParam <- sigmoidalFitVector$slopeParam_Estimate
      incrementTime <- sigmoidalFitVector$incrementTime
      startPoint_x <- sigmoidalFitVector$startPoint_x
      startPoint_y <- sigmoidalFitVector$startPoint_y
      reachMaximum_x <- sigmoidalFitVector$reachMaximum_x
      reachMaximum_y <- sigmoidalFitVector$reachMaximum_y
      h0_est       <- sigmoidalFitVector$h0_Estimate
      if (is.na(fittedXmax)) {
        fittedXmax_sigmoidal <- sigmoidalFitVector$dataScalingParameters.timeRange
      }
      if (!is.na(fittedXmax)) {
        fittedXmax_sigmoidal <- fittedXmax
      }
      if (fittedXmin == 0) {
        fittedXmin_sigmoidal <- 0
      }
      if (fittedXmin != 0) {
        fittedXmin_sigmoidal <- fittedXmin
      }
      time <- seq(fittedXmin_sigmoidal, fittedXmax_sigmoidal,
                  fittedXmax_sigmoidal/1000)
      intensityTheoreticalSigmoidal <- sigmoidalFitFormula_h0(time,
                                                              maximum = maximum_y, slopeParam = slopeParam,
                                                              midPoint = midPoint_x, h0 = h0_est) # appended
      intensityTheoreticalSigmoidalDf <- data.frame(time,
                                                    intensityTheoreticalSigmoidal)
    }
  }
  if (!is.null(doubleSigmoidalFitVector)) {
    if (!doubleSigmoidalFitVector$model == "doublesigmoidal") {
      stop("provided doubleSigmoidalFitVector is not a double sigmoidal fit vector")
    }
    if (!doubleSigmoidalFitVector$isThisaFit) {
      warning("provided double sigmoidal fit vector does not include a fit!")
    }
    if (doubleSigmoidalFitVector$isThisaFit) {
      maximum_x <- doubleSigmoidalFitVector$maximum_x
      maximum_y <- doubleSigmoidalFitVector$maximum_y
      midPoint1_x <- doubleSigmoidalFitVector$midPoint1_x
      midPoint1_y <- doubleSigmoidalFitVector$midPoint1_y
      midPoint1Param <- doubleSigmoidalFitVector$midPoint1Param_Estimate
      midPoint2_x <- doubleSigmoidalFitVector$midPoint2_x
      midPoint2_y <- doubleSigmoidalFitVector$midPoint2_y
      midPointDistanceParam <- doubleSigmoidalFitVector$midPointDistanceParam_Estimate
      slope1 <- doubleSigmoidalFitVector$slope1
      slope1Param <- doubleSigmoidalFitVector$slope1Param_Estimate
      slope2 <- doubleSigmoidalFitVector$slope2
      slope2Param <- doubleSigmoidalFitVector$slope2Param_Estimate
      finalAsymptoteIntensity <- doubleSigmoidalFitVector$finalAsymptoteIntensity
      finalAsymptoteIntensityRatio <- doubleSigmoidalFitVector$finalAsymptoteIntensityRatio_Estimate
      incrementTime <- doubleSigmoidalFitVector$incrementTime
      startPoint_x <- doubleSigmoidalFitVector$startPoint_x
      startPoint_y <- doubleSigmoidalFitVector$startPoint_y
      reachMaximum_x <- doubleSigmoidalFitVector$reachMaximum_x
      reachMaximum_y <- doubleSigmoidalFitVector$reachMaximum_y
      decrementTime <- doubleSigmoidalFitVector$decrementTime
      startDeclinePoint_x <- doubleSigmoidalFitVector$startDeclinePoint_x
      startDeclinePoint_y <- doubleSigmoidalFitVector$startDeclinePoint_y
      endDeclinePoint_x <- doubleSigmoidalFitVector$endDeclinePoint_x
      endDeclinePoint_y <- doubleSigmoidalFitVector$endDeclinePoint_y
      h0_est <- doubleSigmoidalFitVector$h0_Estimate
      if (is.na(fittedXmax)) {
        fittedXmax_doublesigmoidal <- doubleSigmoidalFitVector$dataScalingParameters.timeRange
      }
      if (!is.na(fittedXmax)) {
        fittedXmax_doublesigmoidal = fittedXmax
      }
      if (fittedXmin == 0) {
        fittedXmin_doublesigmoidal <- 0
      }
      if (fittedXmin != 0) {
        fittedXmin_doublesigmoidal <- fittedXmin
      }
      time <- seq(fittedXmin_doublesigmoidal, fittedXmax_doublesigmoidal,
                  fittedXmax_doublesigmoidal/1000)
      intensityTheoreticalDoubleSigmoidal <- doubleSigmoidalFitFormula_h0(time,
                                                                    finalAsymptoteIntensityRatio = finalAsymptoteIntensityRatio,
                                                                    maximum = maximum_y, slope1Param = slope1Param,
                                                                    midPoint1Param = midPoint1Param, slope2Param = slope2Param,
                                                                    midPointDistanceParam = midPointDistanceParam, h0 = h0_est)
      intensityTheoreticalDoubleSigmoidalDf <- data.frame(time,
                                                          intensityTheoreticalDoubleSigmoidal)
    }
  }
  output <- ggplot2::ggplot(dataFrameInput) + ggplot2::geom_point(ggplot2::aes_(x = ~time,
                                                                                y = ~intensity)) + ggplot2::expand_limits(x = 0, y = 0) +
    ggplot2::theme_bw() + ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::xlab(xlabelText) + ggplot2::ylab(ylabelText)
  if (!is.null(sigmoidalFitVector)) {
    if (sigmoidalFitVector$isThisaFit) {
      if (showParameterRelatedLines) {
        if (!sigmoidalFitVector$additionalParameters) {
          stop("to show parameter related lines one needs to run parameterCalculation_h0 for sigmoidalModel ")
        }
        if (sigmoidalFitVector$additionalParameters) {
          output <- output + ggplot2::geom_hline(yintercept = 0,
                                                 colour = "#bdbdbd", size = 0.5, linetype = "longdash") +
            ggplot2::geom_hline(yintercept = maximum_y,
                                colour = "#bdbdbd", size = 0.5, linetype = "longdash") +
            ggplot2::geom_segment(x = startPoint_x,
                                  y = startPoint_y, xend = reachMaximum_x,
                                  yend = reachMaximum_y, colour = "#bdbdbd",
                                  size = 0.5, linetype = "longdash")
        }
      }
      output <- output + ggplot2::geom_point(data = dataFrameInput,
                                             ggplot2::aes_(x = ~time, y = ~intensity)) +
        ggplot2::geom_line(data = intensityTheoreticalSigmoidalDf,
                           ggplot2::aes_(x = ~time, y = ~intensityTheoreticalSigmoidal),
                           color = "orange", size = 1.5)
      if (showParameterRelatedLines) {
        if (sigmoidalFitVector$additionalParameters) {
          output <- output + ggplot2::geom_point(x = midPoint_x,
                                                 y = midPoint_y, colour = "red", size = 6,
                                                 shape = 13)
        }
      }
    }
  }
  if (!is.null(doubleSigmoidalFitVector)) {
    if (doubleSigmoidalFitVector$isThisaFit) {
      if (showParameterRelatedLines) {
        if (!doubleSigmoidalFitVector$additionalParameters) {
          stop("to show parameter related lines one needs to run parameterCalculation_h0 for doubleSigmoidalModel ")
        }
        if (doubleSigmoidalFitVector$additionalParameters) {
          output <- output + ggplot2::geom_hline(yintercept = 0,
                                                 colour = "#bdbdbd", size = 0.5, linetype = "longdash") +
            ggplot2::geom_hline(yintercept = maximum_y,
                                colour = "#bdbdbd", size = 0.5, linetype = "longdash") +
            ggplot2::geom_segment(x = maximum_x, y = finalAsymptoteIntensity,
                                  xend = Inf, yend = finalAsymptoteIntensity,
                                  colour = "#bdbdbd", size = 0.5, linetype = "longdash") +
            ggplot2::geom_segment(x = startPoint_x,
                                  y = startPoint_y, xend = reachMaximum_x,
                                  yend = reachMaximum_y, colour = "#bdbdbd",
                                  size = 0.5, linetype = "longdash") + ggplot2::geom_segment(x = startDeclinePoint_x,
                                                                                             y = startDeclinePoint_y, xend = endDeclinePoint_x,
                                                                                             yend = endDeclinePoint_y, colour = "#bdbdbd",
                                                                                             size = 0.5, linetype = "longdash")
        }
      }
      output <- output + ggplot2::geom_point(data = dataFrameInput,
                                             ggplot2::aes_(x = ~time, y = ~intensity)) +
        ggplot2::geom_line(data = intensityTheoreticalDoubleSigmoidalDf,
                           ggplot2::aes_(x = ~time, y = ~intensityTheoreticalDoubleSigmoidal),
                           color = "orange", size = 1.5)
      if (showParameterRelatedLines) {
        if (doubleSigmoidalFitVector$additionalParameters) {
          output <- output + ggplot2::geom_point(x = maximum_x,
                                                 y = maximum_y, colour = "red", size = 6,
                                                 shape = 13) + ggplot2::geom_point(x = midPoint1_x,
                                                                                   y = midPoint1_y, colour = "red", size = 6,
                                                                                   shape = 13) + ggplot2::geom_point(x = midPoint2_x,
                                                                                                                     y = midPoint2_y, colour = "red", size = 6,
                                                                                                                     shape = 13)
        }
      }
    }
  }
  return(output)
}
