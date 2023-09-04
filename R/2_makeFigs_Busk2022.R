figDir = '../Figs'
library(ErrViewLib)
library(CHNOSZ)
gPars = ErrViewLib::setgPars(type = 'publish')
scalePoints = 0.2

# Busk2022 ####
D = read.table(
  '../Data/BUS2022/qm9_U0_test_Orig.csv',
  sep = ',',
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

S  = D[, "formula"]
R  = D[, "U0"]
C  = D[, "prediction"]
uC = D[, "uncertainty_total"] ^ 0.5
E  = R - C
uE = uC

M = length(uE)

unit = '[eV]'

# subsample
if (M > 25000) {
  sel = sample(1:M, 25000)
  C  = C[sel]
  E  = E[sel]
  uE = uE[sel]
}
Z = E / uE

resVarZ = ErrViewLib::varZCI(Z, method = 'cho')
VarZ    = resVarZ$mean
uVarZ   = resVarZ$sd
print(ErrViewLib::prettyUnc(VarZ, uVarZ, numDig = 1))

# Get mass from formula; use as Input feature
masses = CHNOSZ::mass(S)

# write.csv(
#   cbind(X = masses, V = C, E = E, uE = uE),
#   row.names = FALSE,
#   file = '../BUS2022_Data.csv'
# )

# Figs ####

## Fig 16a ####
png(
  file = file.path(figDir, 'Fig_16a.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotEvsPU(
  uE,
  E,
  logX = TRUE,
  runQuant = TRUE,
  scalePoints = scalePoints,
  ylim = c(-0.3, 0.75),
  xlab = "Uncertainty, uE [eV]",
  ylab = "Error [eV]",
  label = 1,
  gPars = gPars
)
dev.off()

## Fig 16b ####
png(
  file = file.path(figDir, 'Fig_16b.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotEvsPU(
  masses,
  Z,
  type = 'horiz',
  logX = FALSE,
  runQuant = TRUE,
  scalePoints = scalePoints,
  xlim = c(50, 150),
  ylim = c(-5, 5),
  xlab = "Molecular mass [Da]",
  label = 2,
  gPars = gPars
)
dev.off()

## Fig 17a ####
png(
  file = file.path(figDir, 'Fig_17a.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotRelDiag(
  uE,
  E,
  nBin = 10,
  equiPop = TRUE,
  unit = unit,
  logX = TRUE,
  label = 1,
  gPars = gPars
)
dev.off()

## Fig 17b ####
png(
  file = file.path(figDir, 'Fig_17b.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotCalibration(
  E,
  uE,
  score = FALSE,
  label = 2,
  gPars = gPars)
dev.off()
## Fig 17c ####
png(
  file = file.path(figDir, 'Fig_17c.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotCC(
  E,
  uE,
  statS = 'MAE',
  unit = unit,
  label = 3,
  score = FALSE,
  col = 5,
  gPars = gPars
)
dev.off()

## Fig 17d ####
png(
  file = file.path(figDir, 'Fig_17d.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotLZMS(
  uE,
  Z,
  xlab = 'Uncertainty, uE [eV]',
  logX = TRUE,
  nBin = 20,
  equiPop = FALSE,
  score = TRUE,
  label = 4,
  gPars = gPars
)
dev.off()

## Fig 17e ####
png(
  file = file.path(figDir, 'Fig_17e.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotLZMS(
  masses,
  Z,
  xlim = c(50, 150),
  xlab = 'Molecular Mass [Da]',
  logX = FALSE,
  nBin = 20,
  equiPop = FALSE,
  logBin = TRUE,
  score = FALSE,
  label = 5,
  gPars = gPars
)
dev.off()

## Fig 17f ####
png(
  file = file.path(figDir, 'Fig_17f.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotCC(
  E,
  uE,
  ylim = c(0, 0.015),
  showUk = FALSE,
  score = FALSE,
  unit = unit,
  label = 6,
  col = 5,
  gPars = gPars
)
dev.off()

## Fig 18a ####
png(
  file = file.path(figDir, 'Fig_18a.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotCalibration(
  E,
  uE,
  dist = 'Student',
  shape = 4,
  label = 1,
  score = FALSE,
  gPars = gPars
)
dev.off()

## Fig 18b ####
png(
  file = file.path(figDir, 'Fig_18b.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotCC(
  E,
  uE,
  showUk = FALSE,
  statS = 'MAE',
  dist_probref = 'T4',
  unit = unit,
  score = FALSE,
  col = 5,
  label = 2,
  gPars = gPars
)
dev.off()

## Fig 18c ####
png(
  file = file.path(figDir, 'Fig_18c.png'),
  width  = gPars$reso,
  height = gPars$reso
)
ErrViewLib::plotCC(
  E,
  uE,
  ylim = c(0, 0.015),
  dist_probref = 'T4',
  score = FALSE,
  unit = unit,
  showUk = FALSE,
  col = 5,
  label = 3,
  gPars = gPars
)
dev.off()
