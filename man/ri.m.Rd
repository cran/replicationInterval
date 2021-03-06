% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meanFunctions.R
\name{ri.m}
\alias{ri.m}
\title{Replication interval for the mean}
\usage{
ri.m(M, SD = NA, VAR = NA, n, rep.n = NA, prob.level = 0.95)
}
\arguments{
\item{M}{Original study: Mean}

\item{SD}{Original study: Standard deviation. Provide this or variance - not both.}

\item{VAR}{Original study: Variance. Provide this or standard deviation - not both.}

\item{n}{Original study: Sample size}

\item{rep.n}{(optional) Replication study: Sample size. If not specified, n is used.}

\item{prob.level}{(optional 0 to 1 value) Probability level desired (0 to 1). If not specified .95 (i.e., 95 percent) will be used.}
}
\value{
The replication interval and related statistics in list format.
}
\description{
Replication interval for the mean
}
\examples{
ri.m(M=2.53,SD=1.02,n=40,rep.n=80)
}

