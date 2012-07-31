\name{callDenoiseVoxel}
\alias{callDenoiseVoxel}
\title{callDenoiseVoxel}
\description{'callDenoiseVoxel' performs the denoising step of the clustering method.
The function returns a list containing a denoised version of the dataset as an array,
as well as a list for which each element contains a list with the pixel index,
the indexes of its neighboors, the resulting denoised signal, and the variance of the denoised signal   }
\usage{callDenoiseVoxel(fp.data.array, fp.stab.var, fp.alpha = 0.05, 
    fp.proc = "bonferroni")}
\arguments{
  \item{fp.data.array}{an 3D array corresponding to the dataset (the third dimension is the time) }
  \item{fp.stab.var}{a numeric or array indicating the variance of the dataset}
  \item{fp.alpha}{a numeric value indicating the level of the statistical multitest H0}
  \item{fp.proc}{a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
"fdr" method is not implemented yet }
}

\value{returns a list containing:
'details' a list containing for each pixel x time a list of 4 objects:
- 'Vx' a vector containing all the indexes of the denoised pixel's neighboors
- 'Ix' a vector containing the denoised signal
- 'varx' a vector containing the variance of the denoised signal
'denois3D' an array containing the denoised version of the dataset}

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
# #computes the average over time for each pixelxtime before denoising
# avg_before    <- apply(adu340_4small,1:2,mean)
# #computes the average over time for each pixelxtime after denoising
# avg_after     <- apply(denoisres$denois3D,1:2,mean)
# #plotting results
# image(1:nrow(avg_before),1:ncol(avg_before),avg_before)
# x11();image(1:nrow(avg_after),1:ncol(avg_after),avg_after)
}
