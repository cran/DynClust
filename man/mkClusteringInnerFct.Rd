\name{mkClusteringInnerFct}
\alias{mkClusteringInnerFct}
\title{mkClusteringInnerFct}
\description{'mkClusteringInnerFct' returns the clustering function}
\usage{mkClusteringInnerFct(fp.res.denois, fp.var, fp.alpha = 0.05)}
\arguments{
  \item{fp.res.denois}{resulting object obtained after denoising}
  \item{fp.var}{a numeric indicating the variance of the dataset}
  \item{fp.alpha}{a numeric value indicating the level of the statistical multitest H0}
}

\value{returns a list containing:
- 'clusters', the cluster list of two elements:
* Ic a matrix of dimension time x cluster containing the signal in each cluster cluster
* lpix a list containing the indexes of the pixels present in each cluster (as vectors)
- 'lastchange' the index of the cluster last modified in the cluster list
- 'pixtoclust' a vector containing the indexes of the pixels left for clustering    }

\author{Tiffany Lieury,Christophe Pouzat, Yves Rozenholc}






\keyword{internal}

