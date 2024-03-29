% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metapro.R
\name{wFisher}
\alias{wFisher}
\title{wFisher}
\usage{
wFisher(p, weight = NULL, is.onetail = TRUE, eff.sign)
}
\arguments{
\item{p}{A numeric vector of p-values}

\item{weight}{A numeric vector of weight or sample size for each experiment}

\item{is.onetail}{Logical. If set TRUE, p-values are combined without considering the direction of effects, and vice versa. Default: TRUE.}

\item{eff.sign}{A vector of signs of effect sizes. It works when is.onetail = FALSE}
}
\value{
p : Combined p-value

overall.eff.direction : The direction of combined effects.
}
\description{
sample size-weighted Fisher's method
}
\examples{
wFisher(p=c(0.01,0.2,0.8), weight = c(50,60,100),is.onetail=FALSE, eff.sign=c(1,1,1))
}
\references{
Becker BJ (1994). “Combining significance levels.” In Cooper H, Hedges LV (eds.), A handbook of research synthesis, 215–230. Russell Sage, New York.

Fisher RA (1925). Statistical methods for research workers. Oliver and Boyd, Edinburgh.
}
