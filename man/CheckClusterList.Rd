\name{CheckClusterList}
\alias{CheckClusterList}
\title{CheckClusterList}
\description{'CheckClusterList' is used after a cluster is modified or added in the cluster list.
It checks if there exists in the cluster list, clusters which signal is statistically homogeneous
with the last created or modified cluster.
If such clusters exist, the last cluster created or modified is merged with the first homogeneous cluster found.}
\usage{CheckClusterList(fp.clust.list, fp.lastchange, fp.var, fp.alpha = 0.05)}
\arguments{
  \item{fp.clust.list}{a list of clusters containing two elements:
- Ic a matrix of dimension time x cluster containing the signal in each cluster
- lpix a list containing the indexes of the pixels present in each cluster (as vectors) }
  \item{fp.lastchange}{a numeric indicating the index of the cluster last modified in the object 'fp.clust.list'}
  \item{fp.var}{a numeric indicating the variance of the dataset}
  \item{fp.alpha}{a numeric value indicating the level of the statistical multitest H0}
}

\value{returns a list containing:
- 'clusters', the cluster list of two elements:
* Ic a matrix of dimension time x cluster containing the signal in each cluster
* lpix a list containing the indexes of the pixels present in each cluster (as vectors)
- 'lastchange' the index of the cluster last modified in the cluster list }

\author{Tiffany Lieury,Christophe Pouzat, Yves Rozenholc}






\keyword{internal}

