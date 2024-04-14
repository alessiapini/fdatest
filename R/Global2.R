# Two groups comparison with the global test
# arguments:
# data1 and data2: functional data sets of the two groups. They can be either data of class fd or matrices with pointwise evaluations. 
# mu: mean difference under the null hypothesis
# B: number of permutations
# paired: whether the test is paired or not
# dx: optional. If data1 and data2 are fd objects, dx is used as grid step to evaluate them
# stat: optional. Test statistic used for the global test. It can either be: 
#   'Integral' (default) integral of the sample mean difference
#   'Max' max of the sample mean difference
#   'Integral_std' (default) integral of the standardized sample mean difference (t-statistic)
#   'Max_std' (default) integral of the standardized sample mean difference (t-statistic)

Global2 <- function(data1,data2,mu=0,B=1000,paired=FALSE,dx=NULL,stat='Integral'){
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
  J <- dim(coeff1)[2]
  n <- n1+n2
  etichetta_ord <- c(rep(1,n1),rep(2,n2))
  #coeff1 <- coeff1 - matrix(data=mu,nrow=n1,ncol=J)
  
  #print('First step: basis expansion')
  #splines coefficients:
  eval <- coeff <- rbind(coeff1,coeff2)
  p <- dim(coeff)[2]
  
  data.eval <- eval
  #data.eval[1:n1,] <- data.eval[1:n1,] + matrix(data=mu,nrow=n1,ncol=J)
  
  #print('Point-wise tests')
  #univariate permutations
  meandiff2 = (colMeans(coeff[1:n1,]) - colMeans(coeff[(n1+1):n,]))^2
  S1 = cov(coeff[1:n1,])
  S2 = cov(coeff[(n1+1):n,])
  Sp = ((n1-1)*S1 + (n2-1)*S2) / (n1+n2-2)
  T0 <- switch(stat,
               Integral=meandiff2,
               Max=meandiff2,
               Integral_std= meandiff2/diag(Sp),
               Max_std= meandiff2/diag(Sp))
  
  T_coeff <- matrix(ncol=p,nrow=B)
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
    
    meandiff2_perm <- (colMeans(coeff_perm[1:n1,]) - colMeans(coeff_perm[(n1+1):n,]))^2
    S1_perm = cov(coeff_perm[1:n1,])
    S2_perm = cov(coeff_perm[(n1+1):n,])
    Sp_perm = ((n1-1)*S1_perm + (n2-1)*S2_perm) / (n1+n2-2)
    T_coeff[perm,] <- switch(stat,
                             Integral=meandiff2_perm,
                             Max=meandiff2_perm,
                             Integral_std= meandiff2_perm/diag(Sp_perm),
                             Max_std= meandiff2_perm/diag(Sp_perm))
    
  }
  pval <- numeric(p)
  for(i in 1:p){
    pval[i] <- sum(T_coeff[,i]>=T0[i])/B
  }
  
  #combination
  #print('Partition tests')
  all_combs <- rbind(rep(1,p))
  ntests = 1
  adjusted.pval <- numeric(p)
  
  if(stat=='Integral'){
    T0_comb <- sum(T0[which(all_combs[1,]==1)])
    T_comb <- (rowSums(T_coeff[,which(all_combs[1,]==1),drop=FALSE]))
    pval.temp <- mean(T_comb>=T0_comb)
    indexes <- which(all_combs[1,]==1)
    adjusted.pval[indexes] <- pval.temp
  }else if (stat=='Max'){
    T0_comb <- max(T0[which(all_combs[test,]==1)])
    T_comb <- (apply(T_coeff[,which(all_combs[test,]==1)],1,max))
    pval.temp <- mean(T_comb>=T0_comb)
    indexes <- which(all_combs[test,]==1)
    adjusted.pval[indexes] <- pval.temp
  }
  
  
  result = list(
    test = '2pop', mu = mu.eval,
    adjusted_pval = adjusted.pval,
    unadjusted_pval = pval,
    data.eval=data.eval,
    ord_labels = etichetta_ord,
    global_pval = adjusted.pval[1]
  )
  return(result)
}

