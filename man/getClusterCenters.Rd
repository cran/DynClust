\name{getClusterCenters}
\alias{getClusterCenters}
\title{getClusterCenters}
\description{'getClusteredSet' returns the dataset where each pixel x time was replaced by its center}
\usage{getClusterCenters(fp.data.array, fp.res.denois, fp.res.clust)}
\arguments{
  \item{fp.data.array}{an 3D array corresponding to the dataset (the third dimension is the time)}
  \item{fp.res.denois}{a list containing the results of the denoising procedure (cf. 'callDenoiseVoxel')}
  \item{fp.res.clust}{a list containing the results of the clustering procedure (cf. 'ClusteringFct')}
}



\author{Tiffany Lieury,Christophe Pouzat, Yves Rozenholc}





