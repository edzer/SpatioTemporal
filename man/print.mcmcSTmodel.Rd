% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMC.R
\name{print.mcmcSTmodel}
\alias{print.mcmcSTmodel}
\title{Print details for \code{mcmcSTmodel} object}
\usage{
\method{print}{mcmcSTmodel}(x, ...)
}
\arguments{
\item{x}{\code{mcmcSTmodel} object to print information for.}

\item{...}{Ignored additional arguments.}
}
\value{
Nothing
}
\description{
\code{\link[base:print]{print}} method for class \code{mcmcSTmodel}.
}
\examples{
  ##load data
  data(MCMC.mesa.model)
  print(MCMC.mesa.model)

}
\seealso{
Other mcmcSTmodel methods: \code{\link{MCMC.STmodel}},
  \code{\link{density.mcmcSTmodel}},
  \code{\link{plot.density.mcmcSTmodel}},
  \code{\link{plot.mcmcSTmodel}},
  \code{\link{print.summary.mcmcSTmodel}},
  \code{\link{summary.mcmcSTmodel}}
}
\author{
Johan Lindstrom
}
