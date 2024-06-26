% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.fdatest2.R
\name{plot.fdatest2}
\alias{plot.fdatest2}
\title{Plot method for fdatest results on two-population test}
\usage{
\method{plot}{fdatest2}(
  x,
  xrange = c(0, 1),
  alpha1 = 0.05,
  alpha2 = 0.01,
  ylab = "Functional Data",
  main = NULL,
  lwd = 0.5,
  col = c(1, 2),
  ylim = NULL,
  type = "l",
  ...
)
}
\arguments{
\item{x}{The object to be plotted. An object of class "\code{fdatest2}", usually, a result of a call 
to \code{\link{fdatest2}}, \code{\link{IWT2}}, \code{\link{TWT2}}, \code{\link{FDR2}}, or \code{\link{PCT2}}.}

\item{xrange}{Range of the \code{x} axis.}

\item{alpha1}{First level of significance used to select and display significant effects. Default is \code{alpha1 = 0.05}.}

\item{alpha2}{Second level of significance used to select and display significant effects. Default is \code{alpha1 = 0.01}. 
\code{alpha1} and \code{alpha2} are s.t. \code{alpha2 < alpha1}. Otherwise the two values are switched.}

\item{ylab}{Label of \code{y} axis of the plot of functional data. Default is "\code{Functional Data}".}

\item{main}{Plot title.}

\item{lwd}{Line width for the plot of the adjusted p-value function. Default is \code{lwd=1}.}

\item{col}{Colors for the plot of functional data. Default is \code{col = 1}.}

\item{ylim}{Range of the \code{y} axis. Default is \code{NULL}, giving a plot with automatic range for functional data.}

\item{type}{line type for the plot of the adjusted p-value function. Default is type='l'.}

\item{...}{Additional plotting arguments that can be used with function \code{plot}, 
such as \code{\link{graphical parameters}} (see \code{\link{par}}).}
}
\value{
No value returned. 
The function produces a graphical output of the fdatest results:  the plot of the functional data and the one of the adjusted p-values. 
The portions of the domain selected as significant by the test at level \code{alpha1} and \code{alpha2} are highlighted in the plot of the adjusted p-value function and in the one of functional data by gray areas (light and dark gray, respectively).
}
\description{
\code{plot} method for class "\code{fdatest2}".
Plotting function creating a graphical output of the testing procedures
for the mean comparison of two groups: functional data and adjusted p-values are plotted.
}
\examples{
# Importing the NASA temperatures data set
data(NASAtemp)

# Performing the TWT for two populations
TWT.result <- fdatest2(NASAtemp$paris,NASAtemp$milan,method="TWT")

# Plotting the results of the TWT
plot(TWT.result,xrange=c(0,12),main='TWT results for testing mean differences')


# Selecting the significant components at 5\% level
which(TWT.result$adjusted_pval < 0.05)

# Performing the IWT for two populations
IWT.result <- fdatest2(NASAtemp$paris,NASAtemp$milan,method="IWT")

# Plotting the results of the IWT
plot(IWT.result,xrange=c(0,12),main='IWT results for testing mean differences')


# Selecting the significant components at 5\% level
which(IWT.result$adjusted_pval < 0.05)

}
\references{
Pini, A., & Vantini, S. (2017). Interval-wise testing for functional data. \emph{Journal of Nonparametric Statistics}, 29(2), 407-424

Pini, A., Vantini, S., Colosimo, B. M., & Grasso, M. (2018). Domain‐selective functional analysis of variance for supervised statistical profile monitoring of signal data. \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)} 67(1), 55-81.

Abramowicz, K., Hager, C. K., Pini, A., Schelin, L., Sjostedt de Luna, S., & Vantini, S. (2018).
Nonparametric inference for functional‐on‐scalar linear models applied to knee kinematic hop data after injury of the anterior cruciate ligament. \emph{Scandinavian Journal of Statistics} 45(4), 1036-1061.
}
\seealso{
\code{\link{IWTimage}} for the plot of p-values heatmaps (for IWT).
}
