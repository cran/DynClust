\name{MultiTestH0}
\alias{MultiTestH0}
\title{Statistical test of zero mean for dynamics}
\usage{
MultiTestH0(proj.matrix, data.var, thrs = NULL)
}
\arguments{
  \item{proj.matrix}{a matrix whose colummns are tested to
  have zero mean or not.}

  \item{data.var}{a numeric providing the known variance}

  \item{thrs}{a numeric vector of thresholds specified as
  the \eqn{1-\alpha} quantiles of the multiple test under
  the null on each partition. If thrs=NULL (default) return
  the global p-value which is the minimum of all p-values
  obtained on each partition.}
}
\value{
If thrs is provide returns a Boolean vector with length the
number of columns of proj.matrix. Element j is TRUE if the
null hypothesis (no difference with the null vector) is
accepted for column j of proj.matrix. Otherwise, return one
p-value per column.
}
\description{
MultiTestH0 tests if each column vectors of a matrix seen
as a noisy dynamic is of zero mean (\eqn{H_0}) or not. The
multiple statistical test assumes known variance and is
based on a multiple \eqn{\chi^2} test.
}
\author{
Tiffany Lieury, Christophe Pouzat, Yves Rozenholc
}
\references{
Baraud Y., Huet S., Laurent B. \emph{Ann. Stat.} (2003)
Durot C., Rozenholc Y. \emph{Methods Math. Stat.} (2006)
}

