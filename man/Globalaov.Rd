% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Globalaov.R
\name{Globalaov}
\alias{Globalaov}
\title{Global testing procedure for testing functional analysis of variance}
\usage{
Globalaov(
  formula,
  B = 1000,
  method = "residuals",
  dx = NULL,
  recycle = TRUE,
  stat = "Integral"
)
}
\arguments{
\item{formula}{An object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
The output variable of the formula can be either a matrix of dimension \code{c(n,J)} collecting the pointwise evaluations of \code{n} functional data on the same grid of \code{J} points, or a \code{fd} object from the package \code{fda}.}

\item{B}{The number of iterations of the MC algorithm to evaluate the p-values of the permutation tests. The defualt is \code{B=1000}.}

\item{method}{Permutation method used to calculate the p-value of permutation tests. Choose "\code{residuals}" for the permutations of residuals under the reduced model, according to the Freedman and Lane scheme, and "\code{responses}" for the permutation of the responses, according to the Manly scheme.}

\item{dx}{Used only if a \code{fd} object is provided. In this case, \code{dx} is the size of the discretization step of the grid  used to evaluate functional data.
If set to \code{NULL}, a grid of size 100 is used. Default is \code{NULL}.}

\item{stat}{Type of test statistic used for the global test. Possible values are: \code{'Integral'} (default)
for the integral over the domain of the F-test statistic; \code{'Max'} for max over the domain of the F-test statistic.}
}
\value{
\code{Globalaov} returns an object of \code{\link{class}} "\code{IWTaov}". The function \code{summary} is used to obtain and print a summary of the results.
This object is a list containing the following components:
\item{call}{The matched call.}
\item{design_matrix}{The design matrix of the functional-on-scalar linear model.}
\item{unadjusted_pval_F}{Evaluation on a grid of the unadjusted p-value function of the functional F-test.}
\item{adjusted_pval_F}{Evaluation on a grid of the adjusted p-value function of the functional F-test.}
\item{unadjusted_pval_factors}{Evaluation on a grid of the unadjusted p-value function of the functional F-tests on each factor of the analysis of variance (rows).}
\item{adjusted.pval.factors}{adjusted p-values of the functional F-tests on each factor of the analysis of variance (rows) and each basis coefficient (columns).}
\item{Global_pval_F}{Global p-value of the overall test F.}
\item{Global_pval_factors}{Global p-value of test F involving each factor separately.}
\item{data.eval}{Evaluation on a fine uniform grid of the functional data obtained through the basis expansion.}
\item{coeff.regr.eval}{Evaluation on a fine uniform grid of the functional regression coefficients.}
\item{fitted.eval}{Evaluation on a fine uniform grid of the fitted values of the functional regression.}
\item{residuals.eval}{Evaluation on a fine uniform grid of the residuals of the functional regression.}
\item{R2.eval}{Evaluation on a fine uniform grid of the functional R-squared of the regression.}
}
\description{
The function implements the Global Testing procedure for testing mean differences between several
functional populations in a one-way or multi-way functional analysis of variance framework.
Functional data are tested globally and unadjusted and adjusted p-value
functions are provided. The unadjusted p-value function controls the point-wise error rate.
The adjusted p-value function controls the
family-wise error rate weakly. Since this is a global test, the adjusted p-value function is constant.
}
\examples{
# Importing the NASA temperatures data set
data(NASAtemp)
temperature <- rbind(NASAtemp$milan,NASAtemp$paris)
groups <- c(rep(0,22),rep(1,22))

# Performing the test
Global.result <- Globalaov(temperature ~ groups,B=1000)

# Summary of the test results
summary(Global.result)

# Plot of the results
layout(1)
plot(Global.result)

# All graphics on the same device
layout(matrix(1:4,nrow=2,byrow=FALSE))
plot(Global.result,main='NASA data', plot.adjpval = TRUE,xlab='Day',xrange=c(1,365))

}
\references{
Abramowicz, K., Pini, A., Schelin, L., Stamm, A., & Vantini, S. (2022).
“Domain selection and familywise error rate for functional data: A unified framework. 
\emph{Biometrics} 79(2), 1119-1132.

D. Freedman and D. Lane (1983). A Nonstochastic Interpretation of Reported Significance Levels. \emph{Journal of Business & Economic Statistics} 1.4, 292-298.

B. F. J. Manly (2006). Randomization, \emph{Bootstrap and Monte Carlo Methods in Biology}. Vol. 70. CRC Press.
}
\seealso{
See \code{\link{summary.IWTaov}} for summaries and \code{\link{plot.IWTaov}} for plotting the results.
See \code{\link{ITPaovbspline}} for a functional analysis of variance test based on B-spline basis expansion.
See also \code{\link{IWTlm}} to fit and test a functional-on-scalar linear model applying the IWT, and \code{\link{IWT1}}, \code{\link{IWT2}}  for one-population and two-population tests.
}
