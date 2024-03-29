#' @title One population Interval Wise Testing procedure
#'
#' @description The function implements the Interval Wise Testing procedure for testing the center of symmetry of a
#' functional population. Functional data are tested locally and  unadjusted  and adjusted p-value
#' functions are provided. The unadjusted p-value function controls the point-wise error rate. The adjusted p-value function controls the
#' interval-wise error rate.
#'
#' @param data Either pointwise evaluations of the functional data set on a uniform grid, or a \code{fd} object from the package \code{fda}.
#' If pointwise evaluations are provided, \code{data} is a matrix of dimensions \code{c(n,J)}, with \code{J} evaluations on columns and \code{n} units on rows.
#'
#' @param mu The center of symmetry under the null hypothesis. Three possibilities are available for \code{mu}:
#' a constant (in this case, a constant function is used);
#' a \code{J}-dimensional vector containing the evaluations on the same grid which \code{data} are evaluated;
#' a \code{fd} object from the package \code{fda} containing one function.
#' The default is \code{mu=0}.
#'
#' @param B The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=1000}.
#'
#' @param dx Used only if a \code{fd} object is provided. In this case, \code{dx} is the size of the discretization step of the grid  used to evaluate functional data.
#' If set to \code{NULL}, a grid of size 100 is used. Default is \code{NULL}.
#'
#' @param recycle Flag used to decide whether the recycled version of the IWT should be used (see Pini and Vantini, 2017 for details). Default is \code{TRUE}.
#'
#'
#' @return \code{IWT1} returns an object of \code{\link{class}} "\code{IWT1}".
#' An object of class "\code{IWT1}" is a list containing at least the following components:
#' \item{test}{String vector indicating the type of test performed. In this case equal to \code{"1pop"}.}
#' \item{mu}{Evaluation on a grid of the center of symmetry under the null hypothesis (as entered by the user).}
#' \item{unadjusted_pval}{Evaluation on a grid of the unadjusted p-value function.}
#' \item{pval_matrix}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate tests. The element \code{(i,j)} of matrix \code{pval.matrix} contains the p-value of the joint NPC test of the components \code{(j,j+1,...,j+(p-i))}.}
#' \item{adjusted_pval}{Evaluation on a grid of the adjusted p-value function.}
#' \item{data.eval}{Evaluation on a grid of the functional data.}
#'
#' @seealso See also \code{\link{plot.IWT1}} and \code{\link{IWTimage}} for plotting the results.
#'
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#'
#' # Performing the IWT for one population
#' IWT.result <- IWT1(NASAtemp$paris,mu=4)
#'
#' # Plotting the results of the IWT
#' plot(IWT.result,xrange=c(0,12),main='Paris temperatures')
#'
#' # Plotting the p-value heatmap
#' IWTimage(IWT.result,abscissa_range=c(0,12))
#'
#' # Selecting the significant components at 5% level
#' which(IWT.result$adjusted_pval < 0.05)
#'
#' @references A. Pini and S. Vantini (2017).
#' The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. Biometrics 73(3): 835–845.
#'
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data. \emph{Journal of Nonparametric Statistics}, 29(2), 407-424
#'
#' @export

IWT1 <- function(data,mu=0,B=1000,dx=NULL,recycle=TRUE){
  pval.correct <- function(pval.matrix){
    matrice_pval_2_2x <- cbind(pval.matrix,pval.matrix)
    p <- dim(pval.matrix)[2]
    matrice_pval_2_2x <- matrice_pval_2_2x[,(2*p):1]
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
    }
    corrected.pval.matrix <- corrected.pval.matrix[,p:1]
    return(corrected.pval.matrix)
  }

  # data preprocessing
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
    stop("First argument must be either a functional data object or a matrix.")
  }

  if (is.fd(mu)){ # mu is a functional data
    rangeval.mu <- mu$basis$rangeval
    if(sum(rangeval.mu == rangeval)!=2){
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

  n <- dim(coeff)[1]
  p <- dim(coeff)[2]
  data.eval <- coeff <- coeff - matrix(data=mu.eval,nrow=n,ncol=p,byrow=TRUE)

  #univariate permutations
  print('Point-wise tests')
  T0 <- abs(colMeans(coeff))^2  #sample mean
  T_coeff <- matrix(ncol=p,nrow=B)
  for (perm in 1:B){
    signs <- rbinom(n,1,0.5)*2 - 1
    coeff_perm <- coeff*signs
    T_coeff[perm,] <- abs(colMeans(coeff_perm))^2
  }
  pval <- numeric(p)
  for(i in 1:p){
    pval[i] <- sum(T_coeff[,i]>=T0[i])/B
  }

  #combination
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
      print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
    }
  }

  corrected.pval.matrix <- pval.correct(matrice_pval_asymm)
  corrected.pval <- corrected.pval.matrix[1,]

  print('Interval-Wise Testing completed')
  IWT.result <- list(
    test = '1pop', mu = mu.eval,
    adjusted_pval = corrected.pval,
    unadjusted_pval = pval,
    pval_matrix = matrice_pval_asymm,
    data.eval=data.eval)
  class(IWT.result) = 'IWT1'
  return(IWT.result)
}
