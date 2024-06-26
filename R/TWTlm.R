#' @title Threshold-wise testing procedure for testing functional-on-scalar linear models
#'
#' @description The function is used to fit and test functional linear models. 
#' It can be used to carry out regression, and analysis of variance. 
#' It implements the Threshold-wise testing procedure (TWT) for testing the significance of the effects of scalar 
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
#'
#' @return \code{TWTlm} returns an object of \code{\link{class}} "\code{TWTlm}". The function \code{summary} is used to obtain and print a summary of the results.
#' An object of class "\code{ITPlm}" is a list containing at least the following components:
#' \item{call}{call of the function.}
#' \item{design_matrix}{design matrix of the linear model.}
#' \item{unadjusted_pval_F}{unadjusted p-value function of the F test.}
#' \item{adjusted_pval_F}{adjusted p-value function of the F test.}
#' \item{unadjusted_pval_part}{unadjusted p-value functions of the functional t-tests on each covariate, 
#'                         separately (rows) on each domain point (columns).}
#' \item{adjusted_pval_part}{adjusted p-values of the functional t-tests on each covariate (rows) on each domain point (columns).}
#' \item{data.eval}{evaluation of functional data.}
#' \item{coeff.regr.eval}{evaluation of the regression coefficients.}
#' \item{fitted.eval}{evaluation of the fitted values.}
#' \item{residuals.eval}{evaluation of the residuals.}
#' \item{R2.eval}{evaluation of the functional R-suared.}
#'
#' @seealso See \code{\link{summary.TWTlm}} for summaries and \code{\link{plot.TWTlm}} for plotting the results.
#' See also \code{\link{TWTaov}} to fit and test a functional analysis of variance applying the TWT, and  \code{\link{TWT2}} for two-population test.
#'
#'
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#' # Defining the covariates
#' temperature <- rbind(NASAtemp$milan,NASAtemp$paris)
#' groups <- c(rep(0,22),rep(1,22))
#'
#' # Performing the TWT
#' TWT.result <- TWTlm(temperature ~ groups,B=1000)
#' # Summary of the TWT results
#' summary(TWT.result)
#'
#' # Plot of the TWT results
#' layout(1)
#' plot(TWT.result,main='NASA data', plot_adjpval = TRUE,xlab='Day',xrange=c(1,365))
#'
#' # All graphics on the same device
#' layout(matrix(1:6,nrow=3,byrow=FALSE))
#' plot(TWT.result,main='NASA data', plot_adjpval = TRUE,xlab='Day',xrange=c(1,365))
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


TWTlm <- function(formula,
                  B = 1000,
                  method = 'residuals',
                  dx=NULL){
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
  
  print('Threshold-wise tests')
  
  # F-test
  thresholds = c(0,sort(unique(pval_glob)),1)
  adjusted.pval_glob <- pval_glob # we initialize the adjusted p-value as unadjusted one
  pval.tmp <- rep(0,p) # inizialize p-value vector resulting from combined test
  for(test in 1:length(thresholds)){
    #print(paste(test,length(thresholds)))
    # test below threshold
    points.1 = which(pval_glob <= thresholds[test])
    T0_comb = sum(T0_glob[points.1],na.rm=TRUE) # combined test statistic
    T_comb <- (rowSums(T_glob[,points.1,drop=FALSE],na.rm=TRUE))
    pval.test <- mean(T_comb>=T0_comb)
    pval.tmp[points.1] <- pval.test
    # compute maximum
    adjusted.pval_glob = apply(rbind(adjusted.pval_glob,pval.tmp),2,max) 
    
    # test above threshold
    points.2 = which(pval_glob > thresholds[test])
    T0_comb = sum(T0_glob[points.2]) # combined test statistic
    T_comb <- (rowSums(T_glob[,points.2,drop=FALSE],na.rm=TRUE))
    pval.test <- mean(T_comb>=T0_comb)
    pval.tmp[points.2] <- pval.test
    # compute maximum
    adjusted.pval_glob = apply(rbind(adjusted.pval_glob,pval.tmp),2,max) 
    
  }
  
  # F-tests on single factors
  thresholds = c(0,sort(unique(as.numeric(pval_part))),1)
  adjusted.pval_part <- pval_part # we initialize the adjusted p-value as unadjusted one
  
  for(ii in 1:(nvar+1)){
    pval.tmp <- rep(0,p) 
    for(test in 1:length(thresholds)){
      #print(paste(test,length(thresholds)))
      # test below threshold
      points.1 = which(pval_part[ii,] <= thresholds[test])
      T0_comb = sum(T0_part[ii,points.1],na.rm=TRUE) # combined test statistic
      T_comb <- (rowSums(T_part[,ii,points.1,drop=FALSE],na.rm=TRUE))
      pval.test <- mean(T_comb>=T0_comb)
      pval.tmp[points.1] <- pval.test
      # compute maximum
      adjusted.pval_part[ii,] = apply(rbind(adjusted.pval_part[ii,],pval.tmp),2,max) 
      
      # test above threshold
      points.2 = which(pval_part[ii,] > thresholds[test])
      T0_comb = sum(T0_part[ii,points.2]) # combined test statistic
      T_comb <- (rowSums(T_part[,ii,points.2,drop=FALSE],na.rm=TRUE))
      pval.test <- mean(T_comb>=T0_comb)
      pval.tmp[points.2] <- pval.test
      # compute maximum
      adjusted.pval_part[ii,] = apply(rbind(adjusted.pval_part[ii,],pval.tmp),2,max) 
      
    }
    
  }
  
  coeff.regr = regr0$coeff
  coeff.t <- ((coeff.regr))
  
  fitted.regr = regr0$fitted
  fitted.t <- (fitted.regr)
  
  rownames(adjusted.pval_part) = var_names
  rownames(coeff.t) = var_names
  rownames(coeff.regr) = var_names
  rownames(pval_part) = var_names
  
  data.eval <- coeff
  residuals.t = data.eval - fitted.t
  ybar.t = colMeans(data.eval)
  npt <- p
  R2.t = colSums((fitted.t - matrix(data=ybar.t,nrow=n,ncol=npt,byrow=TRUE))^2)/colSums((data.eval - matrix(data=ybar.t,nrow=n,ncol=npt,byrow=TRUE))^2)
  
  print('Threshold-Wise Testing completed')
  
  
  TWT_result <- list(call=cl,
                     design_matrix=design_matrix,
                     unadjusted_pval_F=pval_glob,
                     adjusted_pval_F=adjusted.pval_glob,
                     unadjusted_pval_part=pval_part,
                     adjusted_pval_part=adjusted.pval_part,
                     data.eval=coeff,
                     coeff.regr.eval=coeff.t,
                     fitted.eval=fitted.t,
                     residuals.eval=residuals.t,
                     R2.eval=R2.t)
  class(TWT_result) = 'TWTlm'
  return(TWT_result)
}
