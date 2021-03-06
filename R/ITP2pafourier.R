#' @title Two populations Interval Testing Procedure with Fourier basis (phase-amplitude decomposition)
#'
#' @description The function implements the Interval Testing Procedure for testing the difference between two functional populations evaluated on a uniform grid. Data are represented by means of the Fourier basis expansion with the phase-amplitude decomposition and the significance of the amplitude and phase of each frequency is tested with an interval-wise control of the Family Wise Error Rate.
#'
#' @param data1 Pointwise evaluations of the first population's functional data set on a uniform grid. \code{data1} is a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns and \code{n1} units on rows.
#'
#' @param data2 Pointwise evaluations of the second population's functional data set on a uniform grid. \code{data2} is a matrix of dimensions \code{c(n2,J)}, with \code{J} evaluations on columns and \code{n2} units on rows.
#'
#' @param maxfrequency The maximum frequency to be used in the Fourier basis expansion of data. The default is \code{floor(dim(data1)[2]/2)}, leading to an interpolating expansion.
#'
#' @param B The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=1000}.
#'
#' @param paired A logical indicating whether the test is paired. The default is \code{FALSE}.
#'
#' @return \code{ITP2fourier} returns an object of \code{\link{class}} "\code{ITP2}".
#' An object of class "\code{ITP2}" is a list containing at least the following components:
#' \item{basis}{String vector indicating the basis used for the first phase of the algorithm. Equal to \code{"paFourier"}.}
#' \item{test}{String vector indicating the type of test performed. Equal to \code{"2pop"}.}
#' \item{paired}{Logical indicating whether the test is paired (as entered by the user).}
#' \item{coeff_phase}{Matrix of dimensions \code{c(n,p)} of the \code{p} phases of the Fourier basis expansion. Rows are associated to units and columns to frequencies: the first \code{n1} rows report the coefficients of the first population units and the following \code{n2} rows report the coefficients of the second population units.}
#' \item{coeff_amplitude}{Matrix of dimensions \code{c(n,p)} of the \code{p} amplitudes of the Fourier basis expansion. Rows are associated to units and columns to frequencies: the first \code{n1} rows report the coefficients of the first population units and the following \code{n2} rows report the coefficients of the second population units.}
#' \item{pval_phase}{Unadjusted p-values of the phase tests for each frequency.}
#' \item{pval_amplitude}{Unadjusted p-values of the amplitude tests for each frequency.}
#' \item{pval.matrix_phase}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate tests on phase. The element \code{(i,j)} of matrix \code{pval.matrix_phase} contains the p-value of the joint NPC test of the frequencies \code{(j,j+1,...,j+(p-i))}.}
#' \item{pval.matrix_amplitude}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate tests on amplitude. The element \code{(i,j)} of matrix \code{pval.matrix_amplitude} contains the p-value of the joint NPC test of the frequencies \code{(j,j+1,...,j+(p-i))}.}
#' \item{adjusted.pval_phase}{Adjusted p-values of the phase tests for each frequency.}
#' \item{adjusted.pval_amplitude}{Adjusted p-values of the amplitude tests for each frequency.}
#' \item{labels}{Labels indicating the population membership of each data.}
#' \item{data.eval}{Evaluation on a fine uniform grid of the functional data obtained through the basis expansion.}
#' \item{heatmap.matrix_phase}{Heatmap matrix of p-values for phase (used only for plots).}
#' \item{heatmap.matrix_amplitude}{Heatmap matrix of p-values for amplitude (used only for plots).}
#'
#' @seealso See also \code{\link{plot.ITP2}} and \code{\link{ITPimage}} for plotting the results,
#'  \code{\link{ITP2bspline}} for ITP based on B-spline basis, \code{\link{ITP2fourier}} for Fourier basis
#'  (without phase and amplitude decomposition), and \code{\link{IWT2}} for a two-sample test that is not based on
#'  an a-priori selected basis expansion.
#'
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#'
#' # Performing the ITP
#' ITP.result <- ITP2pafourier(NASAtemp$milan,NASAtemp$paris,maxfrequency=20,B=1000,paired=TRUE)
#'
#' # Plotting the results of the ITP
#' plot(ITP.result,main='NASA data',xrange=c(1,365),xlab='Day')
#'
#' # Plotting the p-value heatmap
#' ITPimage(ITP.result,abscissa.range=c(1,365))
#'
#' # Selecting the significant coefficients
#' which(ITP.result$adjusted.pval < 0.05)
#'
#' @references A. Pini and S. Vantini (2017).
#' The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. Biometrics 73(3): 835–845.
#'
#' @export

ITP2pafourier <-
function(data1,data2,maxfrequency=floor(dim(data1)[2]/2),B=10000,paired=FALSE){
  fisher_cf_L <- function(L){ #fisher on rows of the matrix L
    return(-2*rowSums(log(L)))
  }
  fisher_cf <- function(lambda){ #fisher on vector lambda
    return(-2*sum(log(lambda)))
  }
  sum.dist.sq <- function(m,phi){
    sum((acos(cos(m)*cos(phi)+sin(m)*sin(phi)))^2)
  }
  pval.correct <- function(pval.matrix){
    matrice_pval_2_2x <- cbind(pval.matrix,pval.matrix)
    p <- dim(pval.matrix)[2]
    matrice_pval_2_2x <- matrice_pval_2_2x[,(2*p):1]
    adjusted.pval <- numeric(p)
    for(var in 1:p){
      pval_var <- matrice_pval_2_2x[p,var]
      inizio <- var
      fine <- var #inizio fisso, fine aumenta salendo nelle righe
      for(riga in (p-1):1){
        fine <- fine + 1
        pval_cono <- matrice_pval_2_2x[riga,inizio:fine]
        pval_var <- max(pval_var,pval_cono)
      }
      adjusted.pval[var] <- pval_var
    }
    adjusted.pval <- adjusted.pval[p:1]
    return(adjusted.pval)
  }
  data1 <- as.matrix(data1)
  data2 <- as.matrix(data2)

  n1 <- dim(data1)[1]
  n2 <- dim(data2)[1]
  J <- dim(data1)[2]
  n <- n1+n2
  etichetta_ord <- c(rep(1,n1),rep(2,n2))
  eval <- rbind(data1,data2)

  print('First step: basis expansion')
  #fourier coefficients:
  ak_hat <- NULL
  bk_hat <- NULL
  Period <- J

  for(unit in 1:n){
    #indice <- 1
    data_temp <- eval[unit,]
    abscissa <- 0:(Period-1)
    trasformata <- fft(data_temp)/length(abscissa)
    ak_hat <- rbind(ak_hat,2*Re(trasformata)[1:(maxfrequency+1)])
    bk_hat <- rbind(bk_hat,-2*Im(trasformata)[2:(maxfrequency+1)])
  }


  coeff <- cbind(ak_hat,bk_hat)
  p <- dim(coeff)[2]
  a0 <- coeff[,1]
  ak <- coeff[,2:((p+1)/2)]
  bk <- coeff[,((p+1)/2+1):p]

  #amplitude and phase
  amplitude <- sqrt(ak^2 + bk^2)
  amplitude <- cbind(a0,amplitude)
  phase <- atan2(bk,ak)
  phase <- cbind(pi*(1-sign(a0)),phase)

  #functional data
  K <- p
  if(K %% 2 ==0){
    K <- K+1
  }
  npt <- 1000
  ascissa.smooth <- seq(0, Period, length.out=npt)
  basis <- matrix(0,nrow=npt,ncol=K)
  basis[,1] <- 1/2
  for(i in seq(2,(K-1),2)){
    basis[,i] <- sin(2*pi*(i/2)*ascissa.smooth/Period)
  }
  for(i in seq(3,(K),2)){
    basis[,i] <- cos(2*pi*((i-1)/2)*ascissa.smooth/Period)
  }
  basis.ord <- cbind(basis[,seq(1,K,2)],basis[,seq(2,K-1,2)])
  data.eval <- coeff %*% t(basis.ord)

  #univariate permutations
  print('Second step: joint univariate tests')
  T0_phase <- rep(0,(p+1)/2)
  for(i in 1:((p+1)/2)){
    m1 <- optimize(sum.dist.sq,interval=c(-pi,pi),phi=phase[1:n1,i])$minimum
    m2 <- optimize(sum.dist.sq,interval=c(-pi,pi),phi=phase[(n1+1):n,i])$minimum
    T0_phase[i] <- min(abs(m1-m2),abs(2*pi-m1+m2))
  }
  T0_amplitude <- abs(log(apply(amplitude[1:n1,],2,prod)^(1/n1) / apply(amplitude[(n1+1):n,],2,prod)^(1/n2)))
  T_amplitude <- T_phase <- matrix(nrow=B,ncol=(p+1)/2)
  for (perm in 1:B){
    if(paired==TRUE){
      if.perm <- rbinom(n1,1,0.5)
      coeff_perm <- coeff
      for(couple in 1:n1){
        if(if.perm[couple]==1){
          coeff_perm[c(couple,n1+couple),] <- coeff[c(n1+couple,couple),]
        }
      }
    }else if(paired==FALSE){
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni,]
    }
    ak_perm <- coeff_perm[,2:((p+1)/2)]
    bk_perm <- coeff_perm[,((p+1)/2+1):p]
    a0_perm <- coeff_perm[,1]
    amplitude_perm <- sqrt(ak_perm^2 + bk_perm^2)
    amplitude_perm <- cbind(abs(a0_perm),amplitude_perm)
    T_amplitude[perm,] <- abs(log(apply(amplitude_perm[1:n1,],2,prod)^(1/n1) / apply(amplitude_perm[(n1+1):n,],2,prod)^(1/n2)))

    phase_perm <- atan2(bk_perm,ak_perm)
    phase_perm <- cbind(pi*(1-sign(a0_perm)),phase_perm)
    T_phase_temp <- rep(0,(p+1)/2)
    for(i in 1:((p+1)/2)){
      m1 <- optimize(sum.dist.sq,interval=c(-pi,pi),phi=phase_perm[1:n1,i])$minimum
      m2 <- optimize(sum.dist.sq,interval=c(-pi,pi),phi=phase_perm[(n1+1):n,i])$minimum
      T_phase_temp[i] <- min(abs(m1-m2),abs(2*pi-m1+m2))
    }
    T_phase[perm,] <- T_phase_temp
  }
  pval_phase <- pval_amplitude <- numeric((p+1)/2)
  for(i in 1:(p+1)/2){
    pval_phase[i] <- sum(T_phase[,i]>=T0_phase[i])/B
    pval_amplitude[i] <- sum(T_amplitude[,i]>=T0_amplitude[i])/B
  }


  #combination
  print('Third step: interval-wise combination and correction')
  dim <- (p+1)/2
  q <- numeric(B)
  L_phase <- L_amplitude <- matrix(nrow=B,ncol=dim)
  for(j in 1:dim){
    ordine <- sort.int(T_phase[,j],index.return=T)$ix
    q[ordine] <- (B:1)/(B)
    L_phase[,j] <- q
    ordine <- sort.int(T_amplitude[,j],index.return=T)$ix
    q[ordine] <- (B:1)/(B)
    L_amplitude[,j] <- q
  }

  #asymmetric combination matrix:
  matrice_pval_asymm_phase <- matrix(nrow=dim,ncol=dim)
  matrice_pval_asymm_phase[dim,] <- pval_phase[1:(dim)]
  pval_2x_phase <- c(pval_phase,pval_phase)
  L_2x_phase <- cbind(L_phase,L_phase)
  matrice_pval_asymm_amplitude <- matrix(nrow=dim,ncol=dim)
  matrice_pval_asymm_amplitude[dim,] <- pval_amplitude[1:(dim)]
  pval_2x_amplitude <- c(pval_amplitude,pval_amplitude)
  L_2x_amplitude <- cbind(L_amplitude,L_amplitude)
  for(i in (dim-1):1){
    for(j in 1:dim){
      inf <- j
      sup <- (dim-i)+j
      T0_temp <- fisher_cf(pval_2x_phase[inf:sup])
      T_temp <- fisher_cf_L(L_2x_phase[,inf:sup])
      pval_temp <- sum(T_temp>=T0_temp)/B
      matrice_pval_asymm_phase[i,j] <- pval_temp

      T0_temp <- fisher_cf(pval_2x_amplitude[inf:sup])
      T_temp <- fisher_cf_L(L_2x_amplitude[,inf:sup])
      pval_temp <- sum(T_temp>=T0_temp)/B
      matrice_pval_asymm_amplitude[i,j] <- pval_temp
    }
    print(paste('creating the p-value matrix: end of row ',as.character(dim-i+1),' out of ',as.character(dim),sep=''))
  }

  #symmetric combination matrix
  matrice_pval_symm_phase <- matrice_pval_symm_amplitude <- matrix(nrow=dim,ncol=4*dim)
  for(i in 0:(dim-1)){
    for(j in 1:(2*dim)){
      matrice_pval_symm_phase[dim-i,j+i+dim] <- matrice_pval_asymm_phase[dim-i,(j+1)%/%2]
      matrice_pval_symm_amplitude[dim-i,j+i+dim] <- matrice_pval_asymm_amplitude[dim-i,(j+1)%/%2]
      if(j+i>2*dim-i){
        matrice_pval_symm_phase[dim-i,j+i-dim] <- matrice_pval_asymm_phase[dim-i,(j+1)%/%2]
        matrice_pval_symm_amplitude[dim-i,j+i-dim] <- matrice_pval_asymm_amplitude[dim-i,(j+1)%/%2]
      }
    }
  }

  adjusted.pval_phase <- pval.correct(matrice_pval_asymm_phase)
  adjusted.pval_amplitude <- pval.correct(matrice_pval_asymm_amplitude)

  print('Interval Testing Procedure completed')
  ITP.result <- list(basis='paFourier',test='2pop',paired=as.character(paired),coeff_phase=phase,pval_phase=pval_phase,pval.matrix_phase=matrice_pval_asymm_phase,adjusted.pval_phase=adjusted.pval_phase,
                     coeff_amplitude=amplitude,pval_amplitude=pval_amplitude,pval.matrix_amplitude=matrice_pval_asymm_amplitude,adjusted.pval_amplitude=adjusted.pval_amplitude,labels=etichetta_ord,
                     data.eval=data.eval,heatmap.matrix_phase=matrice_pval_symm_phase,heatmap.matrix_amplitude=matrice_pval_symm_amplitude)
  class(ITP.result) = 'ITP2'
  return(ITP.result)
}
