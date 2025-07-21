test_that("Sigmoidal Fit Function", {
  times <- c(0, 500, 1000)
  sig <- sigmoidalFitFormula_h0(
    x          = times,
    maximum    = 750,
    slopeParam = 0.04,
    midPoint   = 500,
    h0         = 0
  )

  expect_equal(sig[1], 0,             tolerance = 1e-2)  # start at h0 ≈ 0
  expect_equal(sig[2], 750 / 2,       tolerance = 1e-2)  # at t=500: half‐max ≈ 375
  expect_equal(sig[3], 750,           tolerance = 1e-3)  # far out: asymptote ≈ 750
})

test_that("Double-Sigmoidal Fit Function (long tail)", {
  times <- seq(0, 1000, by = 100)

  dsig <- doubleSigmoidalFitFormula_h0(
    x                            = times,
    finalAsymptoteIntensityRatio = 0.5,
    maximum                      = 480,
    slope1Param                  = 0.04,
    midPoint1Param               = 160,
    slope2Param                  = 0.03,
    midPointDistanceParam        = 80,
    h0                           = 0
  )

  expect_equal(max(dsig), 480, tolerance = 1e-2) # reaches max
  expect_equal(dsig[length(dsig)] / max(dsig), 0.5, tolerance = 1e-2) # far out: asymptote ratio ≈ 0.5
})
