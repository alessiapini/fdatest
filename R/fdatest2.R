#' @title Two population local testing procedures
#'
#' @description The function implements local testing procedures for testing mean differences between two
#' functional populations. Functional data are tested locally and unadjusted and adjusted p-value
#' functions are provided. The unadjusted p-value function controls the point-wise error rate. The adjusted p-value function can be computed according
#' to the following methods: 
#' - global testing (controlling the FWER weakly)
#' - interval-wise testing (controlling the interval-wise error rate)
#' - threshold-wise testing (controlling the FWER asymptotically)
#' - partition closed testing (controlling the FWER on a partition)
#' - functional Benjamini Hochberg (controlling the FDR)
#'
#' @param data1 First population's data. Either pointwise evaluations of the functional data set on a uniform grid, or a \code{fd} object from the package \code{fda}.
#' If pointwise evaluations are provided, \code{data2} is a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns and \code{n1} units on rows.
#'
#' @param data2 Second population's data. Either pointwise evaluations of the functional data set on a uniform grid, or a \code{fd} object from the package \code{fda}.
#' If pointwise evaluations are provided, \code{data2} is a matrix of dimensions \code{c(n1,J)}, with \code{J} evaluations on columns and \code{n2} units on rows.
#'
#' @param method A character string specifying the method chosen to adjust the p-value. Should be one of the following: "\code{Global}", "\code{IWT}", "\code{TWT}", "\code{PCT}", "\code{FDR}".
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
#' @param alternative A character string specifying the alternative hypothesis, must be one of "\code{two.sided}" (default), "\code{greater}" or "\code{less}".
#'
#' @param recycle Used only if \code{method}="\code{IWT}". Flag used to decide whether the recycled version of the IWT should be used. Default is \code{TRUE}.
#' 
#' @param partition Used only if \code{method}="\code{PCT}". The partition to be used for PCT procedure. Default is \code{NULL}.
#'
#' @param verbose Used only if \code{method}="\code{IWT}". Logical: if \code{FALSE}, reduces the amount of output. Default is \code{TRUE}.
#'
#' @return \code{fdatest2} returns an object of \code{\link{class}} "\code{fdatest2}" containing at least the following components:
#' \item{test}{String vector indicating the type of test performed. In this case equal to \code{"2pop"}.}
#' \item{mu}{Evaluation on a grid of the functional mean difference under the null hypothesis (as entered by the user).}
#' \item{unadjusted_pval}{Evaluation on a grid of the unadjusted p-value function.}
#' \item{adjusted_pval}{Evaluation on a grid of the adjusted p-value function.}
#' \item{data.eval}{Evaluation on a grid of the functional data.}
#' \item{ord_labels}{Vector of labels indicating the group membership of data.eval}
#'
#' @seealso See also \code{\link{plot.fdatest2}} for plotting the results.
#'
#' @examples
#' # Importing the NASA temperatures data set
#' data(NASAtemp)
#'
#' # Performing the TWT for two populations
#' TWT.result <- fdatest2(NASAtemp$paris,NASAtemp$milan,"TWT")
#'
#' # Plotting the results of the TWT
#' plot(TWT.result,xrange=c(0,12),main='TWT results for testing mean differences')
#'
#'
#' # Selecting the significant components at 5% level
#' which(TWT.result$adjusted_pval < 0.05)
#' 
#' # Performing the IWT for two populations
#' IWT.result <- fdatest2(NASAtemp$paris,NASAtemp$milan,"IWT")
#'
#' # Plotting the results of the IWT
#' plot(IWT.result,xrange=c(0,12),main='IWT results for testing mean differences')
#'
#'
#' # Selecting the significant components at 5% level
#' which(IWT.result$adjusted_pval < 0.05)
#'
#' @references
#' Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
#' “Domain selection and familywise error rate for functional data: A unified framework. 
#' \emph{Biometrics} 79(2), 1119-1132.
#'
#' Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data. \emph{Journal of Nonparametric Statistics}, 29(2), 407-424
#'
#' A. Pini and S. Vantini (2017).
#' The Interval Testing Procedure: Inference for Functional Data Controlling the Family Wise Error Rate on Intervals. Biometrics 73(3): 835–845.
#'
#' Lundtorp Olsen, N., Pini, A., & Vantini, S. (2021).
#' False discovery rate for functional data 
#' \emph{TEST} 30, 784–809.
#'
#' @export

fdatest2 <- function(data1,data2,method,mu=0,B=1000,paired=FALSE,dx=NULL,alternative="two.sided",recycle=TRUE,partition=NULL,verbose=TRUE){
  possible_alternatives <- c("two.sided", "less", "greater")
  if(!(alternative %in% possible_alternatives)){
    stop(paste0('Possible alternatives are ',paste0(possible_alternatives,collapse=', ')))
  }
  
  possible_methods <- c("IWT", "TWT", "PCT", "Global","FDR")
  if(!(method %in% possible_methods)){
    stop(paste0('Available methods are ',paste0(possible_methods,collapse=', ')))
  }
  
  if(method=="PCT" & is.null(partition)){
    stop('PCT method requires to specify a partition of the domain')
  }
  
  result = switch(method,
                  IWT = IWT2(data1=data1,data2=data2,mu=mu,B=B,paired=paired,dx=dx,recycle = recycle,alternative = alternative,verbose=verbose),
                  TWT = TWT2(data1=data1,data2=data2,mu=mu,B=B,paired=paired,dx=dx,alternative = alternative),
                  PCT = PCT2(data1=data1,data2=data2,mu=mu,B=B,paired=paired,dx=dx,alternative = alternative,partition = partition),
                  Global = Global2(data1=data1,data2=data2,mu=mu,B=B,paired=paired,dx=dx,alternative = alternative),
                  FDR = FDR2(data1=data1,data2=data2,mu=mu,B=B,paired=paired,dx=dx,alternative = alternative))
  result$method = method
  return(result)
}

