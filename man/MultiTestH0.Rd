\name{MultiTestH0}
\alias{MultiTestH0}
\title{MultiTestH0}
\description{'MultiTestH0' tests if the vector or matrix 'fp.diff.matrix' given as an argument is significantly different from a null vector or matrix at a level alpha.}
\usage{MultiTestH0(fp.diff.matrix, fp.stab.var, fp.alpha = 0.05, fp.proc = "bonferroni")}
\arguments{
  \item{fp.diff.matrix}{a matrix of signals to be tested with as many columns as signals to test}
  \item{fp.stab.var}{a numeric or array indicating the variance of the dataset}
  \item{fp.alpha}{a numeric value indicating the level of the statistical multitest H0}
  \item{fp.proc}{a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
"fdr" method is not implemented yet}
}

\value{returns TRUE if the test was not significant, FALSE if H1 is true (the mean is centered on 0) at a level alpha
when 'fp.diff.matrix' is a matrix the function returns a vector of boolean values. }

\author{Yves Rozenholc, Christophe Pouzat, Tiffany Lieury}





