figDir = '../Figs'
library(ErrViewLib)
gPars = ErrViewLib::setgPars(type = 'publish')
scalePoints = 0.2

biPlotLZISD = function(M1,
                       M2,
                       set,
                       logX = FALSE,
                       nBin = 10,
                       ylim = NULL,
                       slide = FALSE,
                       method = 'cho',
                       label = 0,
                       xlab = "Uncertainty [eV]",
                       equiPop = TRUE,
                       logBin = FALSE,
                       gPars = ErrViewLib::setgPars()) {
  X1 = M1$uE
  Z1 = M1$E / M1$uE
  intrv1 = ErrViewLib::genIntervals(X1,
                                    nBin,
                                    slide = slide,
                                    equiPop = equiPop,
                                    logBin = logBin)
  xOrd = sort(X1)
  mint1 = c()
  for (i in 1:intrv1$nbr) {
    sel = intrv1$lwindx[i]:intrv1$upindx[i]
    mint1[i] = mean(range(xOrd[sel]))
  }

  X2 = M2$uE
  Z2 = M2$E / M2$uE
  intrv2 = ErrViewLib::genIntervals(X2,
                                    nBin,
                                    slide = slide,
                                    equiPop = equiPop,
                                    logBin = logBin)
  xOrd = sort(X2)
  mint2 = c()
  for (i in 1:intrv2$nbr) {
    sel = intrv2$lwindx[i]:intrv2$upindx[i]
    mint2[i] = mean(range(xOrd[sel]))
  }

  xlim = c(0.9 * min(c(mint1, mint2)), 1.1 * max(c(mint1, mint2)))

  if (is.null(ylim)) {
    out1 = ErrViewLib::plotLZISD(
      X1,
      Z1,
      intrv = intrv1,
      nBin = nBin,
      method = method,
      gPars = gPars,
      plot = FALSE
    )
    out2 = ErrViewLib::plotLZISD(
      X2,
      Z2,
      intrv = intrv2,
      nBin = nBin,
      method = method,
      gPars = gPars,
      plot = FALSE
    )

    ylim = c(min(out1$pcl, out2$pcl), max(out1$pcu, out2$pcu))

  }

  ErrViewLib::plotLZISD(
    X1,
    Z1,
    nBin = nBin,
    intrv = intrv1,
    logX = logX,
    method = method,
    slide = slide,
    xlim = xlim,
    xlab = xlab,
    ylim = ylim,
    col = 2,
    label = label,
    title = set,
    gPars = gPars
  )
  ErrViewLib::plotLZISD(
    X2,
    Z2,
    nBin = nBin,
    intrv = intrv2,
    logX = logX,
    method = method,
    slide = slide,
    col = 5,
    xlim = xlim,
    # Ensures that the average bar lies on the axis
    add = TRUE,
    gPars = gPars
  )
  legend(
    'top',
    bty = 'n',
    legend = c('Uncalibrated', 'Calibrated'),
    col = gPars$cols[c(2, 5)],
    lwd = gPars$lwd,
    lty = 3,
    pch = 16
  )

}
biPlotRelDiag = function(M1,
                         M2,
                         set,
                         logX = FALSE,
                         method = 'cho',
                         equiPop = FALSE,
                         nBin = 10,
                         unit = '',
                         xlim = NULL,
                         logBin = FALSE,
                         label = 0,
                         gPars = ErrViewLib::setgPars()) {
  intrv1 = ErrViewLib::genIntervals(M1$uE,
                                    nBin,
                                    slide = FALSE,
                                    equiPop = equiPop,
                                    logBin = logBin)
  xOrd = sort(M1$uE)
  mint1 = c()
  for (i in 1:intrv1$nbr) {
    sel = intrv1$lwindx[i]:intrv1$upindx[i]
    mint1[i] = sqrt(mean(xOrd[sel] ^ 2))
  }
  intrv2 = ErrViewLib::genIntervals(M2$uE,
                                    nBin,
                                    slide = FALSE,
                                    equiPop = equiPop,
                                    logBin = logBin)
  xOrd = sort(M2$uE)
  mint2 = c()
  for (i in 1:intrv2$nbr) {
    sel = intrv2$lwindx[i]:intrv2$upindx[i]
    mint2[i] = sqrt(mean(xOrd[sel] ^ 2))
  }

  if (is.null(xlim)) {
    # Blank runs for estimation of common plot limits
    out1 = ErrViewLib::plotRelDiag(
      M1$uE,
      M1$E,
      intrv = intrv1,
      method = method,
      logX = logX,
      col = 2,
      gPars = gPars
    )
    ylow = min(out1$pcl)
    yupr = max(out1$pcu)
    out2 = ErrViewLib::plotRelDiag(
      M2$uE,
      M2$E,
      intrv = intrv2,
      method = method,
      logX = logX,
      col = 2,
      gPars = gPars
    )
    ylow = min(ylow, out2$pcl)
    yupr = max(yupr, out2$pcu)
    xlim = range(c(mint1, mint2))
    xlim = ylim = c(min(xlim[1], ylow), max(xlim[2], yupr))

  } else {
    ylim = xlim

  }

  # Plots
  ErrViewLib::plotRelDiag(
    M1$uE,
    M1$E,
    intrv = intrv1,
    method = method,
    nBoot = nBoot,
    logX = logX,
    xlim = xlim,
    ylim = ylim,
    unit = unit,
    col = 2,
    label = label,
    gPars = gPars
  )

  ErrViewLib::plotRelDiag(
    M2$uE,
    M2$E,
    intrv = intrv2,
    method = method,
    nBoot = nBoot,
    xlim = xlim,
    ylim = ylim,
    add = TRUE,
    col = 5,
    gPars = gPars
  )
  legend(
    'topleft',
    bty = 'n',
    legend = c('Uncalibrated', 'Calibrated'),
    col = gPars$cols[c(2, 5)],
    pch = 19,
    lty = 3
  )

}

plotStats = function(df,
                     logX = TRUE,
                     label = 0,
                     title = 'Z-score variance',
                     gPars) {
  for (n in names(gPars))
    assign(n, rlist::list.extract(gPars, n))

  y  = rev(1:nrow(df))
  x  = as.numeric(df[, 2])
  xl = as.numeric(df[, 3])
  xu = as.numeric(df[, 4])
  x1  = as.numeric(df[, 5])
  xl1 = as.numeric(df[, 6])
  xu1 = as.numeric(df[, 7])

  xlim = c(min(c(xl, xl1)), max(c(xu, xu1)))
  if (!logX)
    xlim = c(-1, 1) * max(abs(xlim))

  cex = 1.2 * cex
  par(
    mfrow = c(1, 1),
    mar = c(3, 12, 1, 1),
    mgp = mgp,
    pty = 's',
    tcl = tcl,
    cex = cex,
    lwd = 2 * lwd,
    cex.main = 1
  )

  plot(
    x,
    y,
    type = 'n',
    log = ifelse(logX, 'x', ''),
    pch = '+',
    xlim = xlim,
    xlab = title,
    yaxt = 'n',
    ylab = '',
    col = cols[2]
  )
  mtext(
    df[, 1],
    side = 2,
    at = y,
    las = 1,
    line = 0.2,
    cex = cex
  )
  grid()
  abline(v = c(0, 1),
         col = cols[1],
         lty = 2)
  # points(x1, y, pch = '+', col = cols[5])
  segments(xl1,
           y,
           xu1,
           y,
           lwd = 15 * lwd,
           col = cols[5],
           lend = 2)
  segments(xl,
           y,
           xu,
           y,
           lwd = 10 * lwd,
           col = cols[2],
           lend = 2)
  box()
  legend(
    'bottomright',
    bty = 'n',
    cex = 0.8,
    legend = c('Uncalibrated', 'Calibrated'),
    lty = 1,
    lwd = c(10, 15) * lwd,
    col = cols[c(2, 5)]
  )

  if (label > 0)
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = cex,
      line = 0.3
    )


}

# Start ####

nBin = 20
logX = TRUE
normalize = FALSE
unit = '[eV]'

sets = c('Diffusion', 'Perovskite')[1:2]
methods = c('RF', 'LR', 'GPR_Bayesian')[1:3]

df = data.frame(
  Name = NA,
  muUncal = NA,
  loUncal = NA,
  upUncal = NA,
  muCal   = NA,
  loCal   = NA,
  upCal   = NA
)
df1 = df

fignum = 12
for (method in methods) {
  fignum = fignum + 1
  label = 0
  for (set in sets) {
    fRoot = paste0(set, '_', method)

    case = 'Test_uncal'
    fName = file.path('..', 'Data', 'PAL2022', paste0(fRoot, '_', case, '.csv'))
    if (!file.exists(fName))
      next
    print(fRoot)
    Muncal = read.csv(fName, header = TRUE)
    sel = Muncal$uE > 1e-6 * sd(Muncal$E)
    if (sum(!sel) > 0) {
      print(c(case, ' non-pos unc:', sum(!sel)))
      Muncal = Muncal[sel,]
    }
    print(nrow(Muncal))

    case = 'Test_cal'
    fName = file.path('..', 'Data', 'PAL2022', paste0(fRoot, '_', case, '.csv'))
    if (!file.exists(fName))
      next
    Mcal = read.csv(fName, header = TRUE)
    sel = Mcal$uE > 1e-6 * sd(Mcal$E)
    if (sum(!sel) > 0) {
      print(c(case, ' non-pos unc:', sum(!sel)))
      Mcal = Mcal[sel,]
    }

    # Gather stats
    Z = Muncal$E / Muncal$uE
    zs1 = varZCI(Z, nBoot = length(Z), method = 'cho')
    Z = Mcal$E / Mcal$uE
    zs2 = varZCI(Z, nBoot = length(Z), method = 'cho')
    df = rbind(
      df,
      c(
        Name    = paste0(set, '_', method),
        muUncal = zs1$mean,
        loUncal = zs1$ci[1],
        upUncal = zs1$ci[2],
        muCal   = zs2$mean,
        loCal   = zs2$ci[1],
        upCal   = zs2$ci[2]
      )
    )
    Z = Muncal$E / Muncal$uE
    m1 = mean(Z)
    um1 = sd(Z) / length(Z)
    lo1 = m1 - 1.96 * um1
    up1 = m1 + 1.96 * um1
    Z = Mcal$E / Mcal$uE
    m2 = mean(Z)
    um2 = sd(Z) / length(Z)
    lo2 = m2 - 1.96 * um2
    up2 = m2 + 1.96 * um2
    df1 = rbind(
      df1,
      c(
        Name    = fRoot,
        muUncal = m1,
        loUncal = lo1,
        upUncal = up1,
        muCal   = m2,
        loCal   = lo2,
        upCal   = up2
      )
    )

    ## Figs 13-15 ####
    label = label + 1
    png(
      file = file.path(figDir, paste0('Fig_', fignum, letters[label], '.png')),
      width = gPars$reso,
      height = gPars$reso
    )
    biPlotLZISD(
      Muncal ,
      Mcal,
      fRoot,
      logX = logX,
      nBin = nBin,
      equiPop = FALSE,
      method = 'cho',
      label = label,
      gPars = gPars
    )
    dev.off()

    label = label + 1
    png(
      file = file.path(figDir, paste0('Fig_', fignum, letters[label], '.png')),
      width = gPars$reso,
      height = gPars$reso
    )
    biPlotRelDiag(
      Muncal ,
      Mcal,
      fRoot,
      method = 'cho',
      unit = unit,
      logX = logX,
      nBin = nBin,
      equiPop = FALSE,
      label = label,
      gPars = gPars
    )
    dev.off()

    label = label + 1
    png(
      file = file.path(figDir, paste0('Fig_', fignum, letters[label], '.png')),
      width = gPars$reso,
      height = gPars$reso
    )
    ErrViewLib::plotCC(
      Mcal$E,
      Mcal$uE,
      unit = unit,
      score = FALSE,
      oracle = FALSE,
      showUk = TRUE,
      col = 5,
      label = label,
      gPars = gPars
    )
    dev.off()

  }
}
df  = df[-1,]
df1 = df1[-1,]

## Fig 12a ####
png(
  file = file.path(figDir, 'Fig_12a.png'),
  width  = 1.5 * gPars$reso,
  height = gPars$reso
)
plotStats(
  df1,
  logX = FALSE,
  title = 'Z-score mean',
  label = 1,
  gPars = gPars
)
dev.off()

## Fig 12b ####
png(
  file = file.path(figDir, 'Fig_12b.png'),
  width  = 1.5 * gPars$reso,
  height = gPars$reso
)
plotStats(df,
          title = 'Z-score variance',
          label = 2,
          gPars = gPars)
dev.off()

# Appendix B ####
## Figs 21-22 ####
nBins = c(10, 20, 40)
logX = TRUE

fRoots = c('Diffusion_LR', 'Perovskite_GPR_Bayesian')
fignum = 20

for (fRoot in fRoots) {
  fignum = fignum + 1

  case = 'Test_uncal'
  fName = file.path('../Data/PAL2022', paste0(fRoot, '_', case, '.csv'))
  if (!file.exists(fName))
    next
  print(fRoot)
  Muncal = read.csv(fName, header = TRUE)
  sel = Muncal$uE > 1e-6 * sd(Muncal$E)
  if (sum(!sel) > 0) {
    print(c(case, ' non-pos unc:', sum(!sel)))
    Muncal = Muncal[sel,]
  }
  print(nrow(Muncal))

  case = 'Test_cal'
  fName = file.path('../Data/PAL2022', paste0(fRoot, '_', case, '.csv'))
  if (!file.exists(fName))
    next
  Mcal = read.csv(fName, header = TRUE)
  sel = Mcal$uE > 1e-6 * sd(Mcal$E)
  if (sum(!sel) > 0) {
    print(c(case, ' non-pos unc:', sum(!sel)))
    Mcal = Mcal[sel,]
  }

  label = 0
  for (nBin in nBins) {
    print(nBin)

    title = paste0('n = ', nBin)

    label = label + 1

    png(
      file = file.path(figDir, paste0('Fig_', fignum, letters[label], '.png')),
      width = gPars$reso,
      height = gPars$reso
    )
    biPlotLZISD(
      Muncal ,
      Mcal,
      title,
      equiPop = TRUE,
      label = label,
      ylim = if (fRoot == 'Diffusion_LR')
        c(0.5, 2)
      else
        NULL,
      logX = logX,
      nBin = nBin,
      gPars = gPars
    )
    dev.off()


    png(
      file = file.path(figDir, paste0('Fig_', fignum, letters[label + 3], '.png')),
      width = gPars$reso,
      height = gPars$reso
    )
    biPlotRelDiag(
      Muncal ,
      Mcal,
      title,
      equiPop = TRUE,
      label = label + 3,
      xlim = if (fRoot == 'Diffusion_LR')
        c(0.15, 1.5)
      else
        NULL,
      logX = logX,
      nBin = nBin,
      gPars = gPars,
      unit = unit
    )
    dev.off()

    title = paste0('n = ', nBin, '; adaptive')

    png(
      file = file.path(figDir, paste0('Fig_', fignum, letters[label + 6], '.png')),
      width = gPars$reso,
      height = gPars$reso
    )
    biPlotLZISD(
      Muncal ,
      Mcal,
      title,
      equiPop = FALSE,
      label = label + 6,
      ylim = if (fRoot == 'Diffusion_LR')
        c(0.5, 2)
      else
        NULL,
      logX = logX,
      nBin = nBin,
      gPars = gPars
    )
    dev.off()

    png(
      file = file.path(figDir, paste0('Fig_', fignum, letters[label + 9], '.png')),
      width = gPars$reso,
      height = gPars$reso
    )
    biPlotRelDiag(
      Muncal ,
      Mcal,
      title,
      equiPop = FALSE,
      label = label + 9,
      xlim = if (fRoot == 'Diffusion_LR')
        c(0.15, 1.5)
      else
        NULL,
      logX = logX,
      nBin = nBin,
      gPars = gPars,
      unit = unit
    )
    dev.off()

  }
}
