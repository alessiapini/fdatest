#' @title Two populations Interval Testing Procedure with Fourier basis
#'
#' @description The function implements the Interval Testing Procedure for testing the difference between two functional populations evaluated on a uniform grid. Data are represented by means of the Fourier basis and the significance of each basis coefficient is tested with an interval-wise control of the Family Wise Error Rate.
#'
#' @param data1 Pointwise evaluations of the first population's functional data set on a uniform grid. \code{data1} is a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns and \code{n1} units on rows.
#'
#' @param data2 Pointwise evaluations of the second population's functional data set on a uniform grid. \code{data2} is a matrix of dimensions \code{c(n2,J)}, with \code{J} evaluations on columns and \code{n2} units on rows.
#'
#' @param mu The difference between the first functional population and the second functional population under the null hypothesis. Either a constant (in this case, a constant function is used) or a \code{J}-dimensional vector containing the evaluations on the same grid which \code{data} are evaluated. The default is \code{mu=0}.
#'
#' @param maxfrequency The maximum frequency to be used in the Fourier basis expansion of data. The default is \code{floor(dim(data1)[2]/2)}, leading to an interpolating expansion.
#'
#' @param B The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=1000}.
#'
#' @param paired A logical indicating whether the test is paired. The default is \code{FALSE}.
#'
#' @return \code{ITP2fourier} returns an object of \code{\link{class}} "\code{ITP2}".
#' An object of class "\code{ITP2}" is a list containing at least the following components:
#' \item{basis}{String vector indicating the basis used for the first phase of the algorithm. In this case equal to \code{"Fourier"}.}
#' \item{test}{String vector indicating the type of test performed. In this case equal to \code{"2pop"}.}
#' \item{mu}{Difference between the first functional population and the second functional population under the null hypothesis
#' (as entered by the user).}
#' \item{paired}{Logical indicating whether the test is paired (as entered by the user).}
#' \item{coeff}{Matrix of dimensions \code{c(n,p)} of the \code{p} coefficients of the B-spline basis expansion,
#' with \code{n=n1+n2}. Rows are associated to units and columns to the basis index.
#' The first \code{n1} rows report the coefficients of the first population units and the following \code{n2} rows report the coefficients
#' of the second population units}
#' \item{pval}{Unadjusted p-values for each basis coefficient.}
#' \item{pval.matrix}{Matrix of dimensions \code{c(p,p)} of the p-values of the multivariate tests. The element \code{(i,j)} of matrix \code{pval.matrix} contains the p-value of the joint NPC test of the components \code{(j,j+1,...,j+(p-i))}.}
#' \item{adjusted.pval}{Adjusted p-values for each basis coefficient.}
#' \item{labels}{Labels indicating the population membership of each data.}
#' \item{data.eval}{Evaluation on a fine uniform grid of the functional data obtained through the basis expansion.}
#' \item{heatmap.matrix}{Heatmap matrix of p-values (used only for plots).}
#'
#' @seealso See also \code{\link{plot.ITP2}} and \code{\link{ITPimage}} for plotting the results,
#'  \code{\link{ITP2bspline}} for ITP based on B-spline basis, \code{\link{IWT2}} for a two-sample test that is not based on
#'  an a-priori selected basis expansion.
#'
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#'
#' # Performing the ITP
#' ITP.result <- ITP2fourier(NASAtemp$milan,NASAtemp$paris,maxfrequency=20,B=1000,paired=TRUE)
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
#'
#' @references A. Pini and S. Vantini (2017).
#' The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. Biometrics 73(3): 835–845.
#'
#' @export

ITP2fourier <-
function(data1,data2,mu=0,maxfrequency=floor(dim(data1)[2]/2),B=10000,paired=FALSE){
  fisher_cf_L <- function(L){ #fisher on rows of the matrix L
    return(-2*rowSums(log(L)))
  }
  fisher_cf <- function(lambda){ #fisher on vector lambda
    return(-2*sum(log(lambda)))
  }
  calcola_hotelling_2pop <- function(pop1,pop2,delta.0){ #calcola T^2 pooled per due pop
    mean1 <- colMeans(pop1)
    mean2 <- colMeans(pop2)
    cov1 <- cov(pop1)
    cov2 <- cov(pop2)
    n1 <- dim(pop1)[1]
    n2 <- dim(pop2)[1]
    p <- dim(pop1)[2]
    Sp      <- ((n1-1)*cov1 + (n2-1)*cov2)/(n1+n2-2)
    Spinv   <- solve(Sp)
    T2 <- n1*n2/(n1+n2) * (mean1-mean2-delta.0) %*% Spinv %*% (mean1-mean2-delta.0)
    return(T2)
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
  data2 <- data2 - matrix(data=mu,nrow=n2,ncol=J)

  print('First step: basis expansion')
  #splines coefficients:
  eval <- rbind(data1,data2)
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
  dim <- (p+1)/2

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
  data.eval[(n1+1):n,] <- data.eval[(n1+1):n,] + matrix(data=mu,nrow=n2,ncol=npt)

  print('Second step: joint univariate tests')
  #univariate permutations
  T0 <- numeric(dim)
  for(freq in 2:dim){
    T0[freq] <- calcola_hotelling_2pop(cbind(ak[1:n1,freq-1],bk[1:n1,freq-1]),cbind(ak[(n1+1):n,freq-1],bk[(n1+1):n,freq-1]),c(0,0))
  }
  T0[1] <- abs(mean(coeff[1:n1,1]) - mean(coeff[(n1+1):n,1]))

  T_hotelling <- matrix(nrow=B,ncol=dim)
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
    T_hotelling_temp <- numeric(dim)
    for(freq in 2:dim){
      T_hotelling_temp[freq] <- calcola_hotelling_2pop(cbind(ak_perm[1:n1,freq-1],bk_perm[1:n1,freq-1]),cbind(ak_perm[(n1+1):n,freq-1],bk_perm[(n1+1):n,freq-1]),c(0,0))
    }
    T_hotelling_temp[1] <- abs(mean(coeff_perm[1:n1,1]) - mean(coeff_perm[(n1+1):n,1]))
    T_hotelling[perm,] <- T_hotelling_temp

  }
  pval <- numeric(dim)
  for(i in 1:dim){
    pval[i] <- sum(T_hotelling[,i]>=T0[i])/B
  }

  #combination
  print('Third step: interval-wise combination and correction')
  q <- numeric(B)
  L <- matrix(nrow=B,ncol=dim)
  for(j in 1:dim){
    ordine <- sort.int(T_hotelling[,j],index.return=T)$ix
    q[ordine] <- (B:1)/(B)
    L[,j] <- q
  }

  #asymmetric combination matrix:
  matrice_pval_asymm <- matrix(nrow=dim,ncol=dim)
  matrice_pval_asymm[dim,] <- pval[1:(dim)]
  pval_2x <- c(pval,pval)
  L_2x <- cbind(L,L)
  for(i in (dim-1):1){
    for(j in 1:dim){
      inf <- j
      sup <- (dim-i)+j
      T0_temp <- fisher_cf(pval_2x[inf:sup])
      T_temp <- fisher_cf_L(L_2x[,inf:sup])
      pval_temp <- sum(T_temp>=T0_temp)/B
      matrice_pval_asymm[i,j] <- pval_temp
    }
    print(paste('creating the p-value matrix: end of row ',as.character(dim-i+1),' out of ',as.character(dim),sep=''))
  }

  #symmetric combination matrix
  matrice_pval_symm <- matrix(nrow=dim,ncol=4*dim)
  for(i in 0:(dim-1)){
    for(j in 1:(2*dim)){
      matrice_pval_symm[dim-i,j+i+dim] <- matrice_pval_asymm[dim-i,(j+1)%/%2]
      if(j+i>2*dim-i){
        matrice_pval_symm[dim-i,j+i-dim] <- matrice_pval_asymm[dim-i,(j+1)%/%2]
      }
    }
  }

  adjusted.pval <- pval.correct(matrice_pval_asymm)
  print('Interval Testing Procedure completed')
  ITP.result <- list(basis='Fourier',test='2pop',mu=mu,paired=as.character(paired),coeff=coeff,pval=pval,pval.matrix=matrice_pval_asymm,adjusted.pval=adjusted.pval,labels=etichetta_ord,data.eval=data.eval,heatmap.matrix=matrice_pval_symm)
  class(ITP.result) = 'ITP2'
  return(ITP.result)
}
