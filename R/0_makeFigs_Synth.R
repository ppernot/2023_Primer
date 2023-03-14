figDir = '../Figs'
library(ErrViewLib)
gPars = ErrViewLib::setgPars(type = 'publish')
scalePoints = 0.2

# Functions ####
genQuantTab = function(uE,
                       prob  = seq(0.005, 0.995, by = 0.005),
                       dist  = c('Normal',
                                 'Student',
                                 'Uniform',
                                 'Laplace',
                                 'Normp',
                                 'T4',
                                 'Normp4'),
                       shape = 4) {
  dist = match.arg(dist)

  # Quantile functions
  Normal = function(p, ...)
    qnorm(p)
  Student = function(p, df = 4)
    qt(p, df = df) / sqrt(df / (df - 2))
  Uniform = function(p, ...)
    qunif(p, -sqrt(3), sqrt(3))
  Laplace = function(pr, ...)
    normalp::qnormp(pr, p = 1) / sqrt(gamma(3) / gamma(1))
  Normp = function(pr, df = 4)
    normalp::qnormp(pr, p = df) / sqrt(df ^ (2 / df) * gamma(3 / df) / gamma(1 /
                                                                               df))
  T4 = function(p, ...)
    Student(p, df = 4)
  Normp4 = function(p, ...)
    Normp(p, df = 4)

  quantu = sapply(
    prob,
    FUN = function(x)
      get(dist)(x, df = shape)
  )
  return(list(quant = uE %o% quantu,
              prob  = prob))

}

# Synthetic data ####
M = 5000
prob = (5:1000) / 1000
set.seed(12345)

Case_A = Case_B = Case_C = Case_D = Case_E = Case_F = list()

# Calibrated, consistent, adaptive
Case_A$X     = sort(rnorm(M))
Case_A$R     = 0.1 + Case_A$X ^ 2
Case_A$uE    = 0.1 * abs(Case_A$R)
Case_A$V     = Case_A$R + Case_A$uE * rnorm(M)
Case_A$E     = Case_A$R - Case_A$V
Case_A$dist  = 'Normal'
Case_A$shape = 0
Case_A$Q     = genQuantTab(
  Case_A$uE,
  prob = prob,
  dist = Case_A$dist,
  shape = Case_A$shape
)

# Calibrated, non-consistent, non-adaptive
Case_B       = Case_A
mue          = sqrt(mean(Case_A$uE ^ 2))
Case_B$uE    = mue * rnorm(M, 1, 0.1)
Case_B$Q     = genQuantTab(
  Case_B$uE,
  prob = prob,
  dist = Case_B$dist,
  shape = Case_B$shape
)

# Calibrated VarE; uncalibrated VarZ, non-consistent, non-adaptive
Case_C       = Case_A
Case_C$uE    = sample(Case_A$uE)
Case_C$Q     = genQuantTab(
  Case_C$uE,
  prob = prob,
  dist = Case_C$dist,
  shape = Case_C$shape
)

# Uncalibrated
Case_D       = Case_A
Case_D$uE    = 2 * Case_A$uE
Case_D$Q     = genQuantTab(
  Case_D$uE,
  prob = prob,
  dist = Case_D$dist,
  shape = Case_D$shape
)

# Calibrated, consistent, adaptive; non-normal
Case_E      = Case_A
Case_E$V    = Case_E$R + Case_E$uE * rt(M, df = 4) / sqrt(2)
Case_E$E    = Case_E$R - Case_E$V
Case_E$Q     = genQuantTab(
  Case_E$uE,
  prob = prob,
  dist = Case_E$dist,
  shape = Case_E$shape
)
Case_E$QT4     = genQuantTab(Case_E$uE,
                             prob = prob,
                             dist = 'T4',
                             shape = 4)

# Homoscedastic: calibrated, consistent, adaptive
Case_F       = Case_A
Case_F$uE    = rep(sqrt(mean(Case_F$uE ^ 2)), length(Case_F$uE))
Case_F$V     = Case_F$R + Case_F$uE * rnorm(M)
Case_F$E     = Case_F$R - Case_F$V
Case_F$Q     = genQuantTab(
  Case_F$uE,
  prob = prob,
  dist = Case_F$dist,
  shape = Case_F$shape
)

cases = c('Case_A', 'Case_B', 'Case_C', 'Case_D', 'Case_E', 'Case_F')

probLCP = c(0.25, 0.5, 0.75, 0.95)

# Figs ####

label = 0

for (case in cases) {
  label = label + 1

  X = get(case)

  pcase = case
  substr(pcase, 5, 5) = ' '
  print(pcase)

  Z = X$E / X$uE

  print(var(X$E) / mean(X$uE ^ 2))
  print(var(Z))

  # Get expanded uncertainties from quantiles
  U = matrix(0, nrow = nrow(X$Q$quant), ncol = length(probLCP))
  for (ip in seq_along(probLCP)) {
    pl = (1 - probLCP[ip]) / 2
    ipl = which.min(abs(pl - prob))
    pu = (1 + probLCP[ip]) / 2
    ipu = which.min(abs(pu - prob))
    U[, ip] = 0.5 * (X$Q$quant[, ipu] - X$Q$quant[, ipl])
  }
  if (case == 'Case_E') {
    UT4 = matrix(0,
                 nrow = nrow(X$QT4$quant),
                 ncol = length(probLCP))
    for (ip in seq_along(probLCP)) {
      pl = (1 - probLCP[ip]) / 2
      ipl = which.min(abs(pl - prob))
      pu = (1 + probLCP[ip]) / 2
      ipu = which.min(abs(pu - prob))
      UT4[, ip] = 0.5 * (X$QT4$quant[, ipu] - X$QT4$quant[, ipl])
    }
  }

  ## Fig 1 ####

  png(
    file = paste0(figDir, '/Fig_01', letters[label], '.png'),
    width = gPars$reso,
    height = gPars$reso
  )
  ErrViewLib::plotCalibrationQ(
    X$E,
    X$Q,
    score = FALSE,
    label = label,
    title = pcase,
    gPars = gPars
  )
  dev.off()

  if (var(X$uE) != 0) {
    ## Fig 2 ####
    png(
      file = paste0(figDir, '/Fig_02', letters[label], '.png'),
      width = gPars$reso,
      height = gPars$reso
    )
    ErrViewLib::plotEvsPU(
      X$uE,
      X$E,
      logX = TRUE,
      runQuant = TRUE,
      label = label,
      xlab = "Uncertainty, uE",
      ylim = c(-2, 2),
      scalePoints = scalePoints,
      title = pcase,
      gPars = gPars
    )
    dev.off()

    ## Fig 3 ####
    png(
      file = paste0(figDir, '/Fig_03', letters[label], '.png'),
      width = gPars$reso,
      height = gPars$reso
    )
    ErrViewLib::plotCalibrationQ(
      X$E,
      X$Q,
      cond = X$uE,
      score = FALSE,
      label = label,
      title = pcase,
      gPars = gPars
    )
    dev.off()
    if (case == 'Case_E') {
      png(
        file = paste0(figDir, '/Fig_03', letters[label + 1], '.png'),
        width = gPars$reso,
        height = gPars$reso
      )
      ErrViewLib::plotCalibrationQ(
        X$E,
        X$QT4,
        cond = X$uE,
        score = FALSE,
        label = label + 1,
        title = expression(list("Case E") ~  ~ italic(D) == italic(t)[nu ==
                                                                        4]),
        gPars = gPars
      )
      dev.off()
    }

    ## Fig 4 ####
    png(
      file = paste0(figDir, '/Fig_04', letters[label], '.png'),
      width = gPars$reso,
      height = gPars$reso
    )
    ErrViewLib::plotRelDiag(
      X$uE,
      X$E,
      logX = TRUE,
      score = FALSE,
      label = label,
      title = pcase,
      gPars = gPars
    )
    dev.off()

    # Fig 5 ####
    png(
      file = paste0(figDir, '/Fig_05', letters[label], '.png'),
      width = gPars$reso,
      height = gPars$reso
    )
    ErrViewLib::plotLZISD(
      X$uE,
      Z,
      logX = TRUE,
      xlab = 'Uncertainty, uE',
      score = FALSE,
      label = label,
      method = 'cho',
      title = pcase,
      gPars = gPars
    )
    dev.off()

    ## Fig 6 ####
    png(
      file = paste0(figDir, '/Fig_06', letters[label], '.png'),
      width = gPars$reso,
      height = gPars$reso
    )
    ErrViewLib::plotLCP(
      X$E,
      U,
      ordX = X$uE,
      logX = TRUE,
      prob = probLCP,
      xlab = 'Uncertainty, uE',
      title = pcase,
      label = label,
      gPars = gPars
    )
    dev.off()
    if (case == 'Case_E') {
      png(
        file = paste0(figDir, '/Fig_06', letters[label + 1], '.png'),
        width = gPars$reso,
        height = gPars$reso
      )
      ErrViewLib::plotLCP(
        X$E,
        UT4,
        ordX = X$uE,
        logX = TRUE,
        prob = probLCP,
        xlab = 'Uncertainty, uE',
        title = expression(list("Case E") ~  ~ italic(D) == italic(t)[nu ==
                                                                        4]),
        label = label + 1,
        gPars = gPars
      )
      dev.off()
    }

    ## Fig 7 ####
    png(
      file = paste0(figDir, '/Fig_07', letters[label], '.png'),
      width = gPars$reso,
      height = gPars$reso
    )
    ErrViewLib::plotCC(
      X$E,
      X$uE,
      score = FALSE,
      label = label,
      legLoc = 'bottomright',
      col = 5,
      title = pcase,
      gPars = gPars
    )
    dev.off()
    if (case == 'Case_E') {
      png(
        file = paste0(figDir, '/Fig_07', letters[label + 1], '.png'),
        width = gPars$reso,
        height = gPars$reso
      )
      ErrViewLib::plotCC(
        X$E,
        X$uE,
        score = FALSE,
        dist_probref = 'T4',
        col = 5,
        label = label + 1,
        legLoc = 'bottomright',
        title = expression(list("Case E") ~  ~ italic(D) == italic(t)[nu ==
                                                                        4]),
        gPars = gPars
      )
      dev.off()
    }

  }

  ## Fig 8 ####
  png(
    file = paste0(figDir, '/Fig_08', letters[label], '.png'),
    width = gPars$reso,
    height = gPars$reso
  )
  ErrViewLib::plotEvsPU(
    X$X,
    Z,
    type = 'horiz',
    logX = FALSE,
    runQuant = TRUE,
    xlab = 'Input Feature, X',
    label = label,
    ylim = if (case == 'Case_C')
      c(-30, 30)
    else
      NULL,
    scalePoints = scalePoints,
    title = pcase,
    gPars = gPars
  )
  dev.off()

  ## Fig 9 ####
  png(
    file = paste0(figDir, '/Fig_09', letters[label], '.png'),
    width = gPars$reso,
    height = gPars$reso
  )
  ErrViewLib::plotCalibrationQ(
    X$E,
    X$Q,
    cond = X$X,
    score = FALSE,
    label = label,
    title = pcase,
    gPars = gPars
  )
  dev.off()

  ## Fig 10 ####
  png(
    file = paste0(figDir, '/Fig_10', letters[label], '.png'),
    width = gPars$reso,
    height = gPars$reso
  )
  ErrViewLib::plotLZISD(
    X$X,
    Z,
    logX = FALSE,
    xlab = 'Input Feature, X',
    score = FALSE,
    label = label,
    method = 'cho',
    title = pcase,
    gPars = gPars
  )
  dev.off()

  ## Fig 11 ####
  png(
    file = paste0(figDir, '/Fig_11', letters[label], '.png'),
    width = gPars$reso,
    height = gPars$reso
  )
  ErrViewLib::plotLCP(
    X$E,
    U,
    ordX = X$X,
    logX = FALSE,
    prob = probLCP,
    xlab = 'Input Feature, X',
    title = pcase,
    label = label,
    gPars = gPars
  )
  dev.off()

  # Appendix A ####

  ## Fig 20 ####
  png(
    file = paste0(figDir, '/Fig_20', letters[label], '.png'),
    width = gPars$reso,
    height = gPars$reso
  )
  ErrViewLib::plotEvsPU(
    X$V,
    Z,
    type = 'horiz',
    logX = all(X$V > 0),
    runQuant = TRUE,
    ylim = if (case == 'Case_C')
      c(-30, 30)
    else
      NULL,
    label = label,
    scalePoints = scalePoints,
    title = pcase,
    gPars = gPars
  )
  dev.off()

  png(
    file = paste0(figDir, '/Fig_20', letters[label + 6], '.png'),
    width = gPars$reso,
    height = gPars$reso
  )
  ErrViewLib::plotLZISD(
    X$V,
    Z,
    logX =  all(X$V > 0),
    label = label + 6,
    score = FALSE,
    xlab = 'Predicted Value, V',
    method = 'cho',
    title = pcase,
    gPars = gPars
  )
  dev.off()

}
