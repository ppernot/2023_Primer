figDir = '../Figs'
library(ErrViewLib)
gPars = ErrViewLib::setgPars(type = 'publish')
scalePoints = 0.2

# Data ####
U = read.table(
  '../Data/HU2022/Figure3e_feature/x_feature_distance.csv',
  sep = ',',
  header = FALSE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
E = read.table(
  '../Data/HU2022/Figure3e_feature/y_residual.csv',
  sep = ',',
  header = FALSE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
qHat68 = 0.030882039366943052
qHat95 = 0.07894857763808902
xlab = 'Feature distance [a.u.]'
title = expression(U[F])
title2 = expression(U[list("F,95")] / 1.96)

E = as.vector(unlist(E))
U = as.vector(unlist(U))

U68 = U * qHat68
U95 = U * qHat95

unit = '[eV/sys]'

# Fig 19a ####
png(
  file = file.path(figDir, 'Fig_19a.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotEvsPU(
  U95,
  E,
  logX = TRUE,
  runQuant = TRUE,
  scalePoints = scalePoints,
  ylim = 0.5 * c(-1, 1),
  xlab = expression(paste(U[list("F,95")], " ", group(
    "[", list("eV/sys"), "]"
  ))),
  ylab = expression(paste("Error ", group(
    "[", list("eV/sys"), "]"
  ))),
  label = 1,
  gPars = gPars
)
dev.off()

# Fig 19b ####
png(
  file = file.path(figDir, 'Fig_19b.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotLCP(
  E,
  U = cbind(U68, U95),
  ordX = U,
  nBin = 20,
  equiPop = FALSE,
  prob = c(0.68, 0.95),
  xlab = xlab,
  ylim = c(0.3, 1),
  logX = TRUE,
  slide = FALSE,
  title = title,
  label = 2,
  gPars = gPars
)
dev.off()

# Fig 19c ####
png(
  file = file.path(figDir, 'Fig_19c.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotCC(
  E,
  U95 / 1.96,
  probref = TRUE,
  title = title2,
  unit = unit,
  score = FALSE,
  col = 5,
  label = 3,
  gPars = gPars
)
dev.off()

U_L = read.table(
  '../Data/HU2022/Figure3f_latent/x_latent_distance.csv',
  sep = ',',
  header = FALSE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
U_L = as.vector(unlist(U_L))

qHat68 = 0.01212766581809073
qHat95 = 0.027450564267796015
U68 = U_L * qHat68
U95 = U_L * qHat95

xlab = 'Feature distance [a.u.]'
title2 = expression(U[list("L,95")] / 1.96)


# Fig 19d ####
png(
  file = file.path(figDir, 'Fig_19d.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotEvsPU(
  U95,
  E,
  ylim = 0.5 * c(-1, 1),
  logX = TRUE,
  runQuant = TRUE,
  scalePoints = scalePoints,
  xlab = expression(paste(U[list("L,95")], " ", group(
    "[", list("eV/sys"), "]"
  ))),
  ylab = expression(paste("Error ", group(
    "[", list("eV/sys"), "]"
  ))),
  label = 4,
  gPars = gPars
)
dev.off()

# Fig 19e ####
png(
  file = file.path(figDir, 'Fig_19e.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotLCP(
  E,
  U = cbind(U68, U95),
  ordX = U,
  prob = c(0.68, 0.95),
  nBin = 20,
  equiPop = FALSE,
  xlab = xlab,
  ylim = c(0.3, 1),
  logX = TRUE,
  slide = FALSE,
  label = 5,
  title = expression(U[L]),
  gPars = gPars
)
dev.off()

# Fig 19f ####
png(
  file = file.path(figDir, 'Fig_19f.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotCC(
  E,
  U95 / 1.96,
  probref = TRUE,
  title = title2,
  unit = unit,
  score = FALSE,
  col = 5,
  label = 6,
  gPars = gPars
)
dev.off()
