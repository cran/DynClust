\name{callDenoiseVoxel}
\alias{callDenoiseVoxel}
\title{callDenoiseVoxel}
\description{'callDenoiseVoxel' performs the denoising step of the clustering method.
The function returns a list containing a denoised version of the dataset as an array,
as well as a list for which each element contains a list with the pixel index,
the indexes of its neighboors, the resulting denoised signal, and the variance of the denoised signal   }
\usage{callDenoiseVoxel(fp.data.array, fp.stab.var, fp.depth = 1, fp.alpha = 0.05, 
    fp.mask.size = NULL, fp.mask.auto = F, fp.mask.step = 50, 
    fp.proc = "bonferroni")}
\arguments{
  \item{fp.data.array}{an 3D array corresponding to the dataset (the third dimension is the time).
The number of observations n must be of the form n=2^d}
  \item{fp.stab.var}{a numeric or array indicating the variance of the dataset}
  \item{fp.depth}{a numeric indicating the depth of a voxel}
  \item{fp.alpha}{a numeric value indicating the level of the statistical multitest H0}
  \item{fp.mask.size}{a vector indicating the size (in the x and y dimension) of a circular region defined around the current pixel used to search for neighbours
per default fp.mask.size is defined as c(dim(fp.data.array)[1],dim(fp.data.array)[2])}
  \item{fp.mask.auto}{a boolean indicating if the size (in the x and y dimension) of a circular region defined around the current pixel used to search for neighbours
should grow up until the maximum amount of neighbours are found or that the x and y limits are reached}
  \item{fp.mask.step}{a numeric value indicating the x and y step for the growth of the searching area (only if fp.mask.auto=T) }
  \item{fp.proc}{a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
"fdr" method is not implemented yet }
}
\details{The denoising procedure is applied on each signal. Given a pixel at location x, 

the difference 3D data-array is constructed by subtracting the signal at location x, Fx,

from all the other signals of the data set. The time-homogeneity is statistically tested

for each difference signal at a level alpha. The locations of the pixels for which the null

hypothesis was not rejected are saved, and will be listed as neighbors of the pixel x.

These neighbors are then used to construct a denoised estimation of Fx by averaging the

closest (in space) time-homogeneous signals selected among all the listed neighbors.

Since it is expected that the closer a neighbor is from the current pixel, the most likely

their signals will be coherent, neighbors are ordered by proximity to the current pixel.

Several estimations of the denoised Fx are constructed by averaging the signals of an

increasing number of neighbors. The size of the neighborhood is increasing geometrically

to optimize the computational costs while still obtaining estimates achieving a good statistical

trade-off between bias and variance in step.

The estimation of the denoised Fx using the largest neighborhood size (hence the highest denoising level),

is then compared to all other estimations. If all the estimations are statistically coherent with each other,

the estimation of the denoised signal located at pixel x using the largest neighborhood size becomes

the final denoised Fx. However, if at least one of the estimations is not statistically coherent with

the one having the largest neighborhood size, the estimation of the denoised signal located at pixel x

using the largest neighborhood size is eliminated from the possible estimations.

The process is repeated until all estimates are statistically close to the estimation

computed with the largest neighborhood size.

The number of clusters could be decreased by tuning the level alpha or the variance

during the denoising procedure.

However, by tuning one of these parameters, users expose themselves to the risk of 

obtaining under or over-smoothed estimated signals in the clusters. 

Details about the denoising method and the statistical test for coherence can be found in the references.}
\value{returns a list containing:
'details' a list containing for each pixel x time a list of 4 objects:
- 'Vx' a vector containing all the indexes of the denoised pixel's neighboors
- 'Ix' a vector containing the denoised signal
- 'varx' a vector containing the variance of the denoised signal
'denois3D' an array containing the denoised version of the dataset}
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
# #computes the average over time for each pixelxtime before denoising
# avg_before    <- apply(adu340_4small,1:2,mean)
# #computes the average over time for each pixelxtime after denoising
# avg_after     <- apply(denoisres$denois3D,1:2,mean)
# #plotting results
# image(1:nrow(avg_before),1:ncol(avg_before),avg_before)
# x11();image(1:nrow(avg_after),1:ncol(avg_after),avg_after)
## End(Not run)
}
