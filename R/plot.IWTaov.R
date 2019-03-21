#' @export

plot.IWTaov <- function(x, xrange = c(0,1), 
                        alpha1 = 0.05, alpha2 = 0.01, 
                        plot_adjpval = FALSE,
                        ylim = NULL, col = 1,
                        ylab = 'Functional Data', 
                        main = NULL, lwd = 1, type='l', ...) {
  if (class(x) != "IWTaov") stop("x should be an object of the class IWTaov")
  if (alpha1 < alpha2) {
    temp <- alpha1
    alpha1 <- alpha2
    alpha2 <- temp
  }
  object <- x
  p <- length(object$unadjusted_pval_F)
  n <- dim(t(object$data.eval))[1]
  xmin <- xrange[1]
  xmax <- xrange[2]
  abscissa_pval = seq(xmin, xmax, len = p)
  devAskNewPage(ask = TRUE)  
  main_f <- paste(main, ': Functional Data and F-test')
  main_f <- sub("^ : +", "", main_f)
  
  if (is.null(ylim)) ylim <- range(object$data.eval)
  matplot(abscissa_pval, t(object$data.eval), type = 'l', col = col, main = main_f, 
          ylab = ylab, ylim = ylim, lwd = lwd, ...)
  difference1 <- which(object$adjusted_pval_F < alpha1)
  if (length(difference1) > 0) {
    for (j in 1:length(difference1)) {
      min_rect <- abscissa_pval[difference1[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray90", 
           density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL,
         border = "black")
  }
  difference2 <- which(object$adjusted_pval_F < alpha2)
  if (length(difference2) > 0) {
    for (j in 1:length(difference2)) {
      min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
      max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
      rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray80", 
           density = -2, border = NA)
    }
    rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], col = NULL,
         border = "black")
  }
  matplot(abscissa_pval, t(object$data.eval), type = 'l', 
          col = col, add = TRUE, lwd = lwd, ...)
  
  nvar <- dim(object$adjusted_pval_factors)[1]
  names_all <- colnames(object$design_matrix)
  interaz <- grep(':', names_all)
  for (var in 1:(dim(object$adjusted_pval_factors)[1])) {
    var_name <- rownames(object$adjusted_pval_factors)[var]
    main_t <- paste(main, ': factor', var_name, sep = ' ')
    main_t <- sub("^ : +", "", main_t)
    if (length(grep(':', var_name)) > 0) { # Plot interaction
      var12 <- strsplit(var_name, ':')
      var1 <- var12[[1]][1]
      var2 <- var12[[1]][2]
      dummy_test1 <- grep(var1, names_all)
      dummy_test2 <- grep(var2, names_all)
      dummy_test <- intersect(dummy_test1, dummy_test2)
      colors <- object$design_matrix[, dummy_test]
      if (length(dim(colors)) > 1) {
        colors <- (apply(colors, 1, paste, collapse = ''))
      }
      colors <- as.factor(colors)
    } else { # Plot of a factor
      dummy_test <- grep(var_name, names_all)
      dummy_test <- setdiff(dummy_test, interaz)
      colors <- object$design_matrix[, dummy_test]
      if (length(dim(colors)) > 1) {
        colors <- (apply(colors, 1, paste, collapse = ''))
      }
      colors <- as.factor(colors)
    }
    matplot(abscissa_pval, t(object$data.eval),type = 'l', col = colors, ylim = ylim,
            lwd = 1, main = main_t, ylab = ylab, ...)
    difference1 <- which(object$adjusted_pval_factors[var,] < alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min_rect <- abscissa_pval[difference1[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray90", 
             density = -2, border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], 
           col = NULL, border = "black")
    }
    difference2 <- which(object$adjusted_pval_factors[var,] < alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray80", 
             density = -2, border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], 
           col = NULL, border = "black")
    }
    matlines(abscissa_pval, t(object$data.eval), type='l', col=colors, ...)
    abline(h = 0, lty = 2, col = 1)
  }
  # Plot adjusted p-values
  if (plot_adjpval == TRUE) {
    main_p <- paste(main, ': Adjusted p-values - F-test')
    main_p <- sub("^ : +", "", main_p)
    plot(abscissa_pval, object$adjusted_pval_F, ylim = c(0, 1),
         main = main_p, ylab = 'p-value',type=type,lwd=2, ...)
    difference1 <- which(object$adjusted_pval_F < alpha1)
    if (length(difference1) > 0) {
      for (j in 1:length(difference1)) {
        min_rect <- abscissa_pval[difference1[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray90", 
             density = -2, border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], 
           col = NULL, border = "black")
    }
    difference2 <- which(object$adjusted_pval_F < alpha2)
    if (length(difference2) > 0) {
      for (j in 1:length(difference2)) {
        min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
        max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
        rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray80", 
             density = -2, border = NA)
      }
      rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], 
           col = NULL, border = "black")
    }
    for (j in 0:10) {
      abline(h = j / 10, col = 'lightgray', lty = "dotted")
    }
    lines(abscissa_pval, object$adjusted_pval_F, lwd=2, type=type,...)
    for (var in 1:(dim(object$adjusted_pval_factors)[1])) {
      var_name <- rownames(object$adjusted_pval_factors)[var]
      main_p <- paste(main, ': Adjusted p-values - factor', var_name)
      main_p <- sub("^ : +", "", main_p)
      plot(abscissa_pval, object$adjusted_pval_factors[var, ],
           ylim = c(0, 1), main = main_p, ylab = 'p-value',lwd=2,type=type, ...)
      difference1 <- which(object$adjusted_pval_factors[var,] < alpha1)
      if (length(difference1) > 0) {
        for (j in 1:length(difference1)) {
          min_rect <- abscissa_pval[difference1[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
          max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
          rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray90", 
               density = -2, border = NA)
        }
        rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], 
             col = NULL, border = "black")
      }
      difference2 <- which(object$adjusted_pval_factors[var,] < alpha2)
      if (length(difference2) > 0) {
        for (j in 1:length(difference2)) {
          min_rect <- abscissa_pval[difference2[j]] - (abscissa_pval[2] - abscissa_pval[1])/2
          max_rect <- min_rect + (abscissa_pval[2] - abscissa_pval[1])
          rect(min_rect, par("usr")[3], max_rect, par("usr")[4], col = "gray80", 
               density = -2, border = NA)
        }
        rect(par("usr")[1], par("usr")[3], par("usr")[2],par("usr")[4], 
             col = NULL, border = "black")
      }
      for (j in 0:10) {
        abline(h = j / 10, col = 'lightgray', lty = "dotted")
      }
      lines(abscissa_pval, object$adjusted_pval_factors[var,],type=type,lwd=2, ...)
    }
  }
  devAskNewPage(ask = FALSE)  
}
