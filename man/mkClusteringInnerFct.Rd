\name{mkClusteringInnerFct}
\alias{mkClusteringInnerFct}
\title{mkClusteringInnerFct}
\description{'mkClusteringInnerFct' returns the clustering function}
\usage{mkClusteringInnerFct(fp.res.denois, fp.alpha = 0.05, fp.proc = "bonferroni")}
\arguments{
  \item{fp.res.denois}{resulting object obtained after denoising}
  \item{fp.alpha}{a numeric value indicating the level of the statistical multitest H0}
  \item{fp.proc}{a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
"fdr" method is not implemented yet}
}

\value{returns a list containing:
- 'clusters', the cluster list of three elements:
* Ic a matrix of dimension time x cluster containing the signal in each cluster
* varc a matrix of dimension time x cluster containing variance of the signal in each cluster
* lpix a list containing the indexes of the pixels present in each cluster (as vectors)
- 'lastchange' the index of the cluster last modified in the cluster list
- 'pixtoclust' a vector containing the indexes of the pixels left for clustering    }

\author{Yves Rozenholc, Christophe Pouzat, Tiffany Lieury}






\keyword{internal}

