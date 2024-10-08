#' @title Two population Interval Wise Testing procedure
#'
#' @description The function implements the Interval Wise Testing procedure for testing mean differences between two
#' functional populations. Functional data are tested locally and unadjusted and adjusted p-value
#' functions are provided. The unadjusted p-value function controls the point-wise error rate. The adjusted p-value function controls the
#' interval-wise error rate.
#'
#' @param data1 First population's data. Either pointwise evaluations of the functional data set on a uniform grid, or a \code{fd} object from the package \code{fda}.
#' If pointwise evaluations are provided, \code{data2} is a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns and \code{n1} units on rows.
#'
#' @param data2 Second population's data. Either pointwise evaluations of the functional data set on a uniform grid, or a \code{fd} object from the package \code{fda}.
#' If pointwise evaluations are provided, \code{data2} is a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns and \code{n2} units on rows.
#'
#' @param mu Functional mean difference under the null hypothesis. Three possibilities are available for \code{mu}:
#' a constant (in this case, a constant function is used);
#' a \code{J}-dimensional vector containing the evaluations on the same grid which \code{data} are evaluated;
#' a \code{fd} object from the package \code{fda} containing one function.
#' The default is \code{mu=0}.
#'
#' @param B The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=1000}.
#'
#' @param paired A logical indicating whether a paired test has to be performed. Default is \code{FALSE}.
#'
#' @param dx Used only if a \code{fd} object is provided. In this case, \code{dx} is the size of the discretization step of the grid  used to evaluate functional data.
#' If set to \code{NULL}, a grid of size 100 is used. Default is \code{NULL}.
#'
#' @param recycle A logical used to decide whether the recycled version of the IWT should be used (see Pini and Vantini, 2017 for details). Default is \code{TRUE}.
#'
#' @param alternative A character string specifying the alternative hypothesis, must be one of "\code{two.sided}" (default), "\code{greater}" or "\code{less}".
#' 
#' @param verbose Logical: if \code{FALSE}, reduces the amount of output. Default is \code{TRUE}.
#'
#' @return \code{IWT2} returns an object of \code{\link{class}} "\code{fdatest2}" containing the following components:
#' \item{test}{String vector indicating the type of test performed. In this case equal to \code{"2pop"}.}
#' \item{mu}{Evaluation on a grid of the functional mean difference under the null hypothesis (as entered by the user).}
#' \item{unadjusted_pval}{Evaluation on a grid of the unadjusted p-value function.}
#' \item{pval_matrix}{Matrix of dimensions \code{c(p,p)} of the p-values of the intervalwise tests. The element \code{(i,j)} of matrix \code{pval.matrix} contains the p-value of the test contains the p-value of the test of interval indexed by \code{(j,j+1,...,j+(p-i))}.}
#' \item{adjusted_pval}{Evaluation on a grid of the adjusted p-value function.}
#' \item{data.eval}{Evaluation on a grid of the functional data.}
#' \item{ord_labels}{Vector of labels indicating the group membership of data.eval}
#'
#' @seealso See also \code{\link{plot.fdatest2}} and \code{\link{IWTimage}} for plotting the results.
#'
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#'
#' # Performing the IWT for two populations
#' IWT.result <- IWT2(NASAtemp$paris,NASAtemp$milan)
#'
#' # Plotting the results of the IWT
#' plot(IWT.result,xrange=c(0,12),main='IWT results for testing mean differences')
#'
#' # Plotting the p-value heatmap
#' IWTimage(IWT.result,abscissa_range=c(0,12))
#'
#' # Selecting the significant components at 5% level
#' which(IWT.result$adjusted_pval < 0.05)
#'
#' @references
#' A. Pini and S. Vantini (2017).
#' The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. Biometrics 73(3): 835–845.
#'
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data. \emph{Journal of Nonparametric Statistics}, 29(2), 407-424
#'
#' @export


IWT2 <- function(data1,data2,mu=0,B=1000,paired=FALSE,dx=NULL,recycle=TRUE,alternative="two.sided",verbose=TRUE){
  pval.correct <- function(pval.matrix){
    matrice_pval_2_2x <- cbind(pval.matrix,pval.matrix)
    p <- dim(pval.matrix)[2]
    matrice_pval_2_2x <- matrice_pval_2_2x[,(2*p):1]
    corrected.pval <- numeric(p)
    corrected.pval.matrix <- matrix(nrow=p,ncol=p)
    corrected.pval.matrix[p,] <- pval.matrix[p,p:1]
    for(var in 1:p){
      pval_var <- matrice_pval_2_2x[p,var]
      inizio <- var
      fine <- var #inizio fisso, fine aumenta salendo nelle righe
      for(riga in (p-1):1){
        fine <- fine + 1
        pval_cono <- matrice_pval_2_2x[riga,inizio:fine]
        pval_var <- max(pval_var,pval_cono,na.rm=TRUE)
        corrected.pval.matrix[riga,var] <- pval_var
      }
      corrected.pval[var] <- pval_var
    }
    corrected.pval <- corrected.pval[p:1]
    corrected.pval.matrix <- corrected.pval.matrix[,p:1]
    return(corrected.pval.matrix)
  }

  possible_alternatives <- c("two.sided", "less", "greater")
  if(!(alternative %in% possible_alternatives)){
    stop(paste0('Possible alternatives are ',paste0(possible_alternatives,collapse=', ')))
  }

  # data preprocessing
  if(is.fd(data1)){ # data1 is a functional data object
    rangeval1 <- data1$basis$rangeval
    rangeval2 <- data2$basis$rangeval
    if(is.null(dx)){
      dx <- (rangeval1[2]-rangeval1[1])*0.01
    }
    if(sum(rangeval1 == rangeval2)!=2){
      stop("rangeval of data1 and data2 must coincide.")
    }
    abscissa <- seq(rangeval1[1],rangeval1[2],by=dx)
    coeff1 <- t(eval.fd(fdobj=data1,evalarg=abscissa))
    coeff2 <- t(eval.fd(fdobj=data2,evalarg=abscissa))

  }else if(is.matrix(data1)){
    coeff1 <- data1
    coeff2 <- data2
  }else{
    stop("First argument must be either a functional data object or a matrix.")
  }

  if (is.fd(mu)){ # mu is a functional data
    rangeval.mu <- mu$basis$rangeval
    if(sum(rangeval.mu == rangeval1)!=2){
      stop("rangeval of mu must be the same as rangeval of data.")
    }
    if(is.null(dx)){
      dx <- (rangeval.mu[2]-rangeval.mu[1])*0.01
    }
    abscissa <- seq(rangeval.mu[1],rangeval.mu[2],by=dx)
    mu.eval <- t(eval.fd(fdobj=mu,evalarg=abscissa))
  }else if(is.vector(mu)){
    mu.eval <- mu
  }else{
    stop("Second argument must be either a functional data object or a numeric vector.")
  }

  n1 <- dim(coeff1)[1]
  n2 <- dim(coeff2)[1]
  p <- dim(coeff1)[2]
  n <- n1+n2
  etichetta_ord <- c(rep(1,n1),rep(2,n2))
  coeff1 <- coeff1 - matrix(data=mu.eval,nrow=n1,ncol=p)

  #print('First step: basis expansion')
  #splines coefficients:
  eval <- coeff <- rbind(coeff1,coeff2)

  data.eval <- eval
  data.eval[1:n1,] <- data.eval[1:n1,] + matrix(data=mu.eval,nrow=n1,ncol=p)

  if(verbose)
    print('Point-wise tests')
  #univariate permutations
  meandiff <- colMeans(coeff[1:n1,,drop=FALSE],na.rm=TRUE) - colMeans(coeff[(n1+1):n,,drop=FALSE],na.rm=TRUE)
  sign.diff <- sign(meandiff)
  sign.diff[which(sign.diff==-1)] <- 0
  T0 <- switch(alternative,
               two.sided =  (meandiff)^2,
               greater   =  (meandiff*sign.diff)^2,
               less      =  (meandiff*(sign.diff-1))^2)

  T_coeff <- matrix(ncol=p,nrow=B)
  for (perm in 1:B){
    if(paired){
      if.perm <- rbinom(n1,1,0.5)
      coeff_perm <- coeff
      for(couple in 1:n1){
        if(if.perm[couple]==1){
          coeff_perm[c(couple,n1+couple),] <- coeff[c(n1+couple,couple),]
        }
      }
    }else{
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni,]
    }
    meandiff <- colMeans(coeff_perm[1:n1,,drop=FALSE],na.rm=TRUE) - colMeans(coeff_perm[(n1+1):n,,drop=FALSE],na.rm=TRUE)
    sign.diff <- sign(meandiff)
    sign.diff[which(sign.diff==-1)] <- 0
    T_coeff[perm,] <- switch(alternative,
                             two.sided =  (meandiff)^2,
                             greater   =  (meandiff*sign.diff)^2,
                             less      =  (meandiff*(sign.diff-1))^2)
  }
  pval <- numeric(p)
  for(i in 1:p){
    pval[i] <- sum(T_coeff[,i]>=T0[i])/B
  }

  #combination
  if(verbose)
    print('Interval-wise tests')

  #asymmetric combination matrix:
  matrice_pval_asymm <- matrix(nrow=p,ncol=p)
  matrice_pval_asymm[p,] <- pval[1:p]
  T0_2x <- c(T0,T0)
  T_coeff_2x <- cbind(T_coeff,T_coeff)

  maxrow <- 1
  # con parametro scale
  #maxrow <- p-scale+1

  if(recycle==TRUE){
    for(i in (p-1):maxrow){ # rows
      for(j in 1:p){ # columns
        inf <- j
        sup <- (p-i)+j
        T0_temp <- sum(T0_2x[inf:sup])
        T_temp <- rowSums(T_coeff_2x[,inf:sup])
        pval_temp <- sum(T_temp>=T0_temp)/B
        matrice_pval_asymm[i,j] <- pval_temp
      }
      if(verbose)
        print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
    }
  }else{ # without recycling
    for(i in (p-1):maxrow){ # rows
      for(j in 1:i){ # columns
        inf <- j
        sup <- (p-i)+j
        T0_temp <- sum(T0_2x[inf:sup])
        T_temp <- rowSums(T_coeff_2x[,inf:sup])
        pval_temp <- sum(T_temp>=T0_temp)/B
        matrice_pval_asymm[i,j] <- pval_temp
      }
      if(verbose)
        print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
    }
  }

  corrected.pval.matrix <- pval.correct(matrice_pval_asymm)
  corrected.pval <- corrected.pval.matrix[1,]

  if(verbose)
    print('Interval-Wise Testing completed')
  IWT.result <- list(
    test = '2pop', mu = mu.eval,
    adjusted_pval = corrected.pval,
    unadjusted_pval = pval,
    pval_matrix = matrice_pval_asymm,
    data.eval=data.eval,
    ord_labels = etichetta_ord)
  class(IWT.result) = 'fdatest2'
  return(IWT.result)
}

