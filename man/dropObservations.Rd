% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CV_aux.R
\name{dropObservations}
\alias{dropObservations}
\title{Drop Observations from a STmodel}
\usage{
dropObservations(STmodel, Ind.cv)
}
\arguments{
\item{STmodel}{Model object from which to drop observations.}

\item{Ind.cv}{A logical vector with one element per observation in \cr
\code{STmodel$obs}. Observations marked with the \code{TRUE}
will be dropped from the data structure. Use \code{\link{createCV}} to
create the logical vector.}
}
\value{
Returns the \code{STmodel} without the observations marked by
  \code{Ind.cv}. Only observed locations are retained.
}
\description{
Drops observations from \code{STmodel}, removing marked observations
along with the corresponding locations and recomputes a number of relevant
elements.
}
\examples{
##load data
data(mesa.model)

##Mark 30\% of observations
I <- runif(dim(mesa.model$obs)[1])<.3
##drop these observations
mesa.model.new <- dropObservations(mesa.model, I)

##This reduces the remaining number of observations
print(mesa.model)
print(mesa.model.new)

\dontshow{
  if( (length(mesa.model$obs$obs)<=length(mesa.model.new$obs$obs)) ||
     (length(mesa.model.new$obs$obs)!=sum(!I)) ){
    stop("drop.observations 1: Observations not dropped")
  }
}
##create cross validation structure
Icv <- createCV(mesa.model, groups=10)

##drop observations from the second CV group
mesa.model.new <- dropObservations(mesa.model, Icv==2)

##This reduces the remaining number of observations (and locations)
print(mesa.model)
print(mesa.model.new)

\dontshow{
  if( (length(mesa.model$obs$obs)<=length(mesa.model.new$obs$obs)) ||
     (length(mesa.model.new$obs$obs)!=sum(Icv!=2)) ||
     (dim(mesa.model$locations)[1] <= dim(mesa.model.new$locations)[1]) ){
    stop("drop.observations 2: Observations not dropped")
  }
}

}
\seealso{
Other STmodel functions: \code{\link{createCV}},
  \code{\link{createDataMatrix}},
  \code{\link{createSTmodel}},
  \code{\link{estimateBetaFields}},
  \code{\link{loglikeSTdim}}, \code{\link{loglikeST}},
  \code{\link{predictNaive}}, \code{\link{processLUR}},
  \code{\link{processLocation}}, \code{\link{updateCovf}},
  \code{\link{updateTrend.STdata}}

Other cross-validation functions: \code{\link{computeLTA}},
  \code{\link{createCV}}, \code{\link{estimateCV.STmodel}},
  \code{\link{predictNaive}}
}
\author{
Johan Lindstrom
}
