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
\details{As for the denoising procedure, the classes are built such that within each class

the clustered signals are statistically coherent. The clustering procedure is a recursive

algorithm, and the main loop of which can be decomposed into the following steps:

- the list of all the pixels, Lv, is ordered by increasing neighborhood sizes;

- a class is defined from a single pixel: at step 1, the class is defined using

the first pixel of Lv. The coherence between the signal in the current pixel and

all other pixels in Lv is statistically tested. The pixels for which the null hypothesis

was not rejected are removed from Lv and included in the current cluster. The estimated

signal of the cluster is computed as the averaged signal over all the pixels of the cluster.

- the estimated signal of the new class is statistically tested for coherence with the

estimated signal of all other existing clusters;

- statistically coherent clusters are merged.

The algorithm stops when Lv is empty and when the clusters cannot be merged anymore.

The number of clusters could be decreased by tuning the level alpha or the variance

during the denoising procedure.

However, by tuning one of these parameters, users expose themselves to the risk of 

obtaining under or over-smoothed estimated signals in the clusters.

Details about the statistical test for coherence can be found in the references.}
\value{returns a list containing:
- 'details' a list containing for each cluster a list of 3 objects:
* Ic a matrix of dimension time x cluster containing the signal in each cluster
* varc a matrix of dimension time x cluster containing variance of the signal in each cluster
* lpix a list containing the indexes of the pixels present in each cluster (as vectors)
- 'clust3D' an array containing the clustered version of the dataset
- 'clustmap' a matrix indicating using their creation index the location of each cluster}
\references{Rozenholc, Y. and Reiss, M. (2012) _Preserving time structures while denoising a dynamical image_,

Mathematical Methods for Signal and Image Analysis and Representation (Chapter 12), 

Florack, L. and Duits, R. and Jongbloed, G. and van~Lieshout, M.-C. and Davies, L.

Ed., Springer-Verlag, Berlin}
\author{Tiffany Lieury,Christophe Pouzat, Yves Rozenholc}




\examples{
## Not run:
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
# denoisres     <- callDenoiseVoxel(adu340_4small,FT_varhat)
# #launches the clustering step on the dataset with a statistical level of 5\%
# clustres      <- ClusteringFct(denoisres)
# x11()
# matplot(clustres$details$Ic,t="l",lty=1,lwd=2,bty="n")
# x11()
# par(mar=rep(0,4))
# image(clustres$clustmap,col=rainbow(max(clustres$clustmap)))
# x11()
# par(mar=rep(0,4))
# image(apply(adu340_4small,1:2,mean),col=grey(seq(0,1,len=2^8)))
## End(Not run)
}
