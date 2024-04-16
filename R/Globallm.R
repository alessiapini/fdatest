#' @title Global testing procedure for testing functional-on-scalar linear models
#'
#' @description The function is used to fit and test functional linear models. 
#' It can be used to carry out regression, and analysis of variance. 
#' It implements the global testing procedure for testing the significance of the effects of scalar 
#' covariates on a functional population. 
#'
#' @param formula An object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted. 
#' Example: y ~ A + B
#'            where: y is a matrix of dimension n * p containing the point-wise evaluations of the n functional data on p points
#'            or an object of class \code{fd} (see \code{fda} package) containing the functional data set
#'            A, B are n-dimensional vectors containing the values of two covariates.
#'            Covariates may be either scalar or factors.
#'
#' @param B The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=1000}.
#'
#' @param method Permutation method used to calculate the p-value of permutation tests. Choose "\code{residuals}" for the permutations of residuals under the reduced model, according to the Freedman and Lane scheme, and "\code{responses}" for the permutation of the responses, according to the Manly scheme.
#'
#' @param dx step size for the point-wise evaluations of functional data. dx is only used ia an object of 
#' class 'fd' is provided as response in the formula.
#' 
#' @param stat Type of test statistic used for the global test. Possible values are: \code{'Integral'} (default)
#' for the integral over the domain of the F- and t-test statistics; \code{'Max'} for max over the domain of the F- and t-test statistics. 

#'
#' @return \code{Globallm} returns an object of \code{\link{class}} "\code{IWTlm}". The function \code{summary} is used to obtain and print a summary of the results.
#' This object is a list containing the following components:
#' \item{call}{call of the function.}
#' \item{design_matrix}{design matrix of the linear model.}
#' \item{unadjusted_pval_F}{unadjusted p-value function of the F test.}
#' \item{adjusted_pval_F}{adjusted p-value function of the F test.}
#' \item{unadjusted_pval_part}{unadjusted p-value functions of the functional t-tests on each covariate, 
#'                         separately (rows) on each domain point (columns).}
#' \item{adjusted_pval_part}{adjusted p-values of the functional t-tests on each covariate (rows) on each domain point (columns).}
#' \item{Global_pval_F}{Global p-value of the overall test F.}
#' \item{Global_pval_part}{Global p-value of t-test involving each covariate separately.}
#' \item{data.eval}{evaluation of functional data.}
#' \item{coeff.regr.eval}{evaluation of the regression coefficients.}
#' \item{fitted.eval}{evaluation of the fitted values.}
#' \item{residuals.eval}{evaluation of the residuals.}
#' \item{R2.eval}{evaluation of the functional R-suared.}
#'
#' @seealso See \code{\link{summary.IWTlm}} for summaries and \code{\link{plot.IWTlm}} for plotting the results.
#' See \code{\link{ITPlmbspline}} for a functional linear model test based on an a-priori selected B-spline basis expansion.
#' See also \code{\link{IWTaov}} to fit and test a functional analysis of variance applying the IWT, and \code{\link{IWT1}}, \code{\link{IWT2}} for one-population and two-population tests.
#'
#'
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#' # Defining the covariates
#' temperature <- rbind(NASAtemp$milan,NASAtemp$paris)
#' groups <- c(rep(0,22),rep(1,22))
#'
#' # Performing the IWT
#' Global.result <- Globallm(temperature ~ groups,B=1000)
#' # Summary of the IWT results
#' summary(Global.result)
#'
#' # Plot of the IWT results
#' layout(1)
#' plot(Global.result,main='NASA data', plot_adjpval = TRUE,xlab='Day',xrange=c(1,365))
#'
#' # All graphics on the same device
#' layout(matrix(1:6,nrow=3,byrow=FALSE))
#' plot(Global.result,main='NASA data', plot_adjpval = TRUE,xlab='Day',xrange=c(1,365))
#'
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified framework. 
#' \emph{Biometrics} 79(2), 1119-1132.
#'
#' D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported Significance Levels. \emph{Journal of Business & Economic Statistics} 1(4), 292-298.
#'
#' B. F. J. Manly (2006). Randomization, \emph{Bootstrap and Monte Carlo Methods in Biology}. Vol. 70. CRC Press.
#'
#' @export


Globallm <- function(formula,
                  B = 1000,
                  method = 'residuals',
                  dx=NULL,
                  recycle=TRUE,
                  stat='Integral'){
  extract_residuals <- function(regr){
    return(regr$residuals)
  }
  extract_fitted <- function(regr){
    return(regr$fitted)
  }
  
  env <- environment(formula)
  variables = all.vars(formula)
  y.name = variables[1]
  covariates.names <- colnames(attr(terms(formula),"factors"))
  #data.all = model.frame(formula)
  cl <- match.call()
  data <- get(y.name,envir = env)
  if(is.fd(data)){ # data is a functional data object
    rangeval <- data$basis$rangeval
    if(is.null(dx)){
      dx <- (rangeval[2]-rangeval[1])*0.01
    }
    abscissa <- seq(rangeval[1],rangeval[2],by=dx)
    coeff <- t(eval.fd(fdobj=data,evalarg=abscissa))
  }else if(is.matrix(data)){
    coeff <- data
  }else{
    stop("First argument of the formula must be either a functional data object or a matrix.")
  }
  
  dummynames.all <- colnames(attr(terms(formula),"factors"))
  formula.const <- deparse(formula[[3]],width.cutoff = 500L) #extracting the part after ~ on formula. this will not work if the formula is longer than 500 char
  
  formula.discrete <- as.formula(paste('coeff ~',formula.const),env=environment())
  design_matrix = model.matrix(formula.discrete)
  mf = model.frame(formula.discrete)
  
  nvar <- dim(design_matrix)[2] - 1
  var_names <- colnames(design_matrix)
  p <- dim(coeff)[2]
  n <- dim(coeff)[1]
  # Univariate permutations
  regr0 <- lm.fit(design_matrix, coeff)
  # Test statistics
  Sigma <- chol2inv(regr0$qr$qr)
  resvar <- colSums(regr0$residuals ^ 2) / regr0$df.residual
  se <- sqrt(matrix(diag(Sigma), nrow = nvar + 1, ncol = p, byrow = FALSE) 
             * matrix(resvar, nrow = nvar + 1, ncol = p, byrow = TRUE))
  T0_part <- abs(regr0$coeff / se)^2
  if (nvar > 0) {
    T0_glob <- colSums((regr0$fitted - matrix(colMeans(regr0$fitted),
                                              nrow = n, ncol = p, 
                                              byrow = TRUE)) ^ 2) / (nvar * resvar)
  } else {
    method <- 'responses' 
    T0_glob <- numeric(p)
    T0_part <- t(as.matrix(T0_part))
  }
  # Compute residuals
  if (method == 'residuals') {
    # n residuals for each coefficient of basis expansion (1:p) 
    # and for each partial test + global test (nvar+1) 
    # Saved in array of dim (nvar+1,n,p)
    # Extracting the part after ~ on formula. 
    # This will not work if the formula 
    # is longer than 500 char
    formula_const <- deparse(formula[[3]], width.cutoff = 500L)
    design_matrix_names2 <- design_matrix
    var_names2 <- var_names
    coeffnames <- paste('coeff[,', as.character(1:p),']', sep = '')
    formula_temp <- coeff ~ design_matrix
    mf_temp <- cbind(model.frame(formula_temp)[-((p + 1):(p + nvar + 1))], 
                     as.data.frame(design_matrix[, -1]))
    if (length(grep('factor', formula_const)) > 0) {
      index_factor <- grep('factor', var_names)
      replace_names <- paste('group', (1:length(index_factor)), sep = '')
      var_names2[index_factor] <- replace_names
      colnames(design_matrix_names2) <- var_names2
    }
    residui <- array(dim=c(nvar + 1, n, p))
    fitted_part <- array(dim = c(nvar + 1, n, p)) 
    formula_coeff_part <- vector('list', nvar + 1)
    regr0_part <- vector('list',nvar + 1)
    # The first one is the intercept. Treated as special case after loop
    for (ii in 2:(nvar + 1)) { 
      var_ii <- var_names2[ii]
      variables_reduced <- var_names2[-c(1, which(var_names2 == var_ii))] 
      if (nvar > 1) {
        formula_temp <- paste(variables_reduced, collapse = ' + ')
      } else {
        # Removing the unique variable -> reduced model only has intercept ter
        formula_temp <- '1' 
      }
      formula_temp2 <- coeff ~ design_matrix_names2
      mf_temp2 <- cbind(model.frame(formula_temp2)[-((p + 1):(p + nvar + 1))], 
                        as.data.frame(design_matrix_names2[,-1]))
      formula_coeff_temp <- paste(coeffnames, '~', formula_temp) 
      formula_coeff_part[[ii]] <- sapply(formula_coeff_temp, as.formula)
      regr0_part[[ii]] <- lapply(formula_coeff_part[[ii]], lm, data = mf_temp2)
      residui[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_residuals))
      fitted_part[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_fitted))
    }
    ii <- 1 # intercept
    formula_temp <- paste(formula_const, ' -1', sep = '')
    formula_coeff_temp <- paste(coeffnames, '~', formula_temp)
    formula_coeff_part[[ii]] <- sapply(formula_coeff_temp, as.formula)
    regr0_part[[ii]] <- lapply(formula_coeff_part[[ii]], lm, data = mf_temp)
    residui[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_residuals))
    fitted_part[ii, , ] <- simplify2array(lapply(regr0_part[[ii]], extract_fitted))
  }
  print('Point-wise tests')
  
  # CMC algorithm
  T_glob <- matrix(ncol = p,nrow = B)
  T_part <- array(dim = c(B, nvar + 1, p))
  for (perm in 1:B) {
    # the F test is the same for both methods
    if (nvar > 0) {
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni, ]
    }else{ # Test on intercept permuting signs
      signs <- rbinom(n, 1, 0.5) * 2 - 1
      coeff_perm <- coeff * signs
    }
    regr_perm <- lm.fit(design_matrix, coeff_perm)
    Sigma <- chol2inv(regr_perm$qr$qr)
    resvar <- colSums(regr_perm$residuals ^ 2) / regr_perm$df.residual
    if (nvar > 0) {
      T_glob[perm, ] <- colSums((regr_perm$fitted - matrix(colMeans(regr_perm$fitted), 
                                                           nrow = n, ncol = p,
                                                           byrow=TRUE)) ^ 2)/ (nvar * resvar)
    }
    # Partial tests: differ depending on the method
    if (method == 'responses') {
      se <- sqrt(matrix(diag(Sigma), nrow = nvar + 1, ncol = p, byrow = FALSE) 
                 * matrix(resvar, nrow = nvar + 1, ncol = p, byrow = TRUE))
      T_part[perm, , ] <- abs(regr0$coeff / se)^2
    } else if (method == 'residuals'){
      residui_perm <- residui[, permutazioni, ]
      regr_perm_part <- vector('list', nvar + 1)
      for (ii in 1:(nvar + 1)) { 
        coeff_perm <- fitted_part[ii, , ] + residui_perm[ii, , ]  
        regr_perm <- lm.fit(design_matrix, coeff_perm)
        Sigma <- chol2inv(regr_perm$qr$qr)
        resvar <- colSums(regr_perm$residuals ^ 2) / regr_perm$df.residual
        se <- sqrt(matrix(diag(Sigma), nrow = nvar + 1 , ncol = p, byrow = FALSE) 
                   * matrix(resvar, nrow = nvar + 1, ncol = p, byrow = TRUE))
        T_part[perm, ii, ] <- abs(regr_perm$coeff / se)[ii, ]^2
      }
    }
  }
  pval_glob <- numeric(p)
  pval_part <- matrix(nrow = nvar + 1, ncol = p)
  for (i in 1:p) {
    pval_glob[i] <- sum(T_glob[, i] >= T0_glob[i]) / B
    pval_part[, i] <- colSums(T_part[, , i] >= matrix(T0_part[, i],nrow = B, ncol = nvar + 1,byrow = TRUE)) / B
  }
  
  print('Global test')
  
  
  if(stat=='Integral'){
    T0_temp <- sum(T0_glob[1:p])
    T_temp <- rowSums(T_glob[,1:p])
    Global_pval_F <- sum(T_temp>=T0_temp)/B
    
    Global_pval_part = numeric(nvar + 1)
    for(ii in 1:(nvar + 1)){
      T0_temp <- sum(T0_part[ii,1:p])
      T_temp <- rowSums(T_part[,ii,1:p])
      Global_pval_part[ii] <- sum(T_temp>=T0_temp)/B
    }
    
    corrected.pval_glob = rep(Global_pval_F,p)
    
    corrected.pval_part = matrix(nrow=nvar+1,ncol=p)
    for(ii in 1:(nvar + 1)){
      corrected.pval_part[ii,] = rep(Global_pval_part[ii],p)
    }
  }else if(stat=='Max'){
    T0_temp <- max(T0_glob)
    T_temp <- apply(T_glob,1,max)
    Global_pval_F <- sum(T_temp>=T0_temp)/B
    
    Global_pval_part = numeric(nvar + 1)
    for(ii in 1:(nvar + 1)){
      T0_temp <- max(T0_part[ii,])
      T_temp <- apply(T_part[,ii,],1,max)
      Global_pval_part[ii] <- sum(T_temp>=T0_temp)/B
    }
    
    corrected.pval_glob = rep(Global_pval_F,p)
    corrected.pval_part = matrix(nrow=nvar+1,ncol=p)
    for(ii in 1:(nvar + 1)){
      corrected.pval_part[ii,] = rep(Global_pval_part[ii],p)
    }
    
  }
  
  
  
  coeff.regr = regr0$coeff
  coeff.t <- ((coeff.regr))
  
  fitted.regr = regr0$fitted
  fitted.t <- (fitted.regr)
  
  rownames(corrected.pval_part) = var_names
  rownames(coeff.t) = var_names
  rownames(coeff.regr) = var_names
  rownames(pval_part) = var_names
  
  data.eval <- coeff
  residuals.t = data.eval - fitted.t
  ybar.t = colMeans(data.eval)
  npt <- p
  R2.t = colSums((fitted.t - matrix(data=ybar.t,nrow=n,ncol=npt,byrow=TRUE))^2)/colSums((data.eval - matrix(data=ybar.t,nrow=n,ncol=npt,byrow=TRUE))^2)
  
  print('Global Testing completed')
  
  
  ITP_result <- list(call=cl,
                     design_matrix=design_matrix,
                     unadjusted_pval_F=pval_glob,
                     adjusted_pval_F=corrected.pval_glob,
                     unadjusted_pval_part=pval_part,
                     adjusted_pval_part=corrected.pval_part,
                     Global_pval_F = Global_pval_F,
                     Global_pval_part = Global_pval_part,
                     data.eval=coeff,
                     coeff.regr.eval=coeff.t,
                     fitted.eval=fitted.t,
                     residuals.eval=residuals.t,
                     R2.eval=R2.t)
  class(ITP_result) = 'IWTlm'
  return(ITP_result)
}
