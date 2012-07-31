\name{ClusteringFct}
\alias{ClusteringFct}
\title{ClusteringFct}
\description{'ClusteringFct' clusters together the denoised pixels x time
with  statistically homogeneous signals at a level 'fp.alpha'}
\usage{ClusteringFct(fp.res.listdenois, fp.alpha = 0.05, fp.proc = "bonferroni")}
\arguments{
  \item{fp.res.listdenois}{a list containing the results of the denoising procedure (cf. 'callDenoiseVoxel')}
  \item{fp.alpha}{a numeric value indicating the level of the statistical multitest H0}
  \item{fp.proc}{a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
"fdr" method is not implemented yet}
}

\value{returns a list containing:
- 'details' a list containing for each cluster a list of 3 objects:
* Ic a matrix of dimension time x cluster containing the signal in each cluster
* varc a matrix of dimension time x cluster containing variance of the signal in each cluster
* lpix a list containing the indexes of the pixels present in each cluster (as vectors)
- 'clust3D' an array containing the clustered version of the dataset
- 'clustmap' a matrix indicating using their creation index the location of each cluster}

\author{Yves Rozenholc, Christophe Pouzat, Tiffany Lieury}




\examples{
# library(DynClust)
# data("adu340_4small",package="DynClust")
# #gain of the CCD camera
# #necessary in order to compute the variance of the dataset
# #estimated on calibration experiments
# G             <- 0.146
# #readout variance
# sro2          <- (16.4)^2
# #dataset's variance
# FT_varhat     <- G*adu340_4small+G^2*sro2
# #launches the denoising step on the dataset with a statistical level of 5\%
# denoisres     <- callDenoiseVoxel(adu340_4small,FT_varhat,.05,"bonferroni")
# #launches the clustering step on the dataset with a statistical level of 5\%
# clustres      <- ClusteringFct(denoisres,.05)
# x11()
# matplot(clustres$details$Ic,t="l",lty=1,lwd=2,bty="n")
# x11()
# par(mar=rep(0,4))
# image(clustres$clustmap,col=rainbow(max(clustres$clustmap)))
# x11()
# par(mar=rep(0,4))
# image(apply(adu340_4small,1:2,mean),col=grey(seq(0,1,len=2^8)))
}
