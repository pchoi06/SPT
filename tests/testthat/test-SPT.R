test_that('SPT with double-sigmoidal input', {
  set.seed(4747)
  reps <- 5
  time_raw <- rep(seq(0, 300, by = 50), each = reps)

  y_true <- doubleSigmoidalFitFormula_spt(
    x                            = time_raw,
    finalAsymptoteIntensityRatio = 0.5,
    maximum                      = 480,
    slope1Param                  = 0.04,
    midPoint1Param               = 160,
    slope2Param                  = 0.03,
    midPointDistanceParam        = 80,
    h0                           = 0
  )

  intensity_raw <- y_true + rnorm(length(y_true), sd = 100)
  dataInput <- data.frame(time = time_raw, intensity = intensity_raw)

  fitObj <- SPT(dataInput,
                             threshold_minimum_for_intensity_maximum = 0.3,
                             threshold_intensity_range = 0.1)

  # check that final decision is correct
  expect_equal(fitObj$decisionProcess$decision, "double_sigmoidal")

  # more tests needed here
})
