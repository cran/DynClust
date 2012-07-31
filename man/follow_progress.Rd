\name{follow_progress}
\alias{follow_progress}
\title{follow progress}
\description{'follow_progress' displays the progression of an iterative process
and returns fp.prog a vector of values ranging from 0 to 100
each value corresponding to a progression bar left for display }
\usage{follow_progress(fp.idx.curr = 0, fp.max.iter = 0, fp.prog = NULL, 
    fp.step = 1)}
\arguments{
  \item{fp.idx.curr}{a numeric indicating which is the current iteration }
  \item{fp.max.iter}{a numeric indicating the total number of iterations}
  \item{fp.prog}{a vector of values ranging between 0 to 100 with a step of fp.step
where each value corresponds to a progression bar to be displayed
per default (when fp.prog=NULL) fp.prog = seq(0,100,by=fp.step)}
  \item{fp.step}{a numeric indicating the step necessary for constructing the fp.prog vector}
}

\value{returns an updated version of fp.prog}

\author{Yves Rozenholc, Christophe Pouzat, Tiffany Lieury}




\examples{
# ratiobars  <-  follow_progress()
#  for(i in 1:5000) ratiobars  <-  follow_progress(i,5000,ratiobars) 
}

\keyword{internal}

