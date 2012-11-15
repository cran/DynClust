\name{mkFCTdenoiseVoxel}
\alias{mkFCTdenoiseVoxel}
\title{mkFCTdenoiseVoxel}
\description{'mkFCTdenoiseVoxel' returns a function used to denoise the pixel x time at index fp.idx.pix,
leaving out of the denoising process the pixels x time at indexes fp.idx.pix.
This function is used to initialise the constant variables
when the function 'callDenoiseVoxel' is called }
\usage{mkFCTdenoiseVoxel(fp.data.coord, fp.data.matrix, fp.stab.var, 
    fp.depth = 1, fp.alpha = 0.05, fp.mask.size = NULL, fp.mask.auto = F, 
    fp.mask.step = 50, fp.proc = "bonferroni")}
\arguments{
  \item{fp.data.coord}{a matrix containing all the pixel x time's (x,y) coordinates
with x in the first column and y in the second }
  \item{fp.data.matrix}{a matrix of pixels x time, each row corresponding to a timepoint and each column to a pixel}
  \item{fp.stab.var}{a numeric indicating the variance of the dataset}
  \item{fp.depth}{a numeric indicating the depth of a voxel}
  \item{fp.alpha}{a numeric value indicating the level of the statistical multitest H0}
  \item{fp.mask.size}{a vector indicating the size (in the x and y dimension) of a circular region defined around the current pixel used to search for neighbours
per default fp.mask.size is defined as c(dim(fp.data.array)[1],dim(fp.data.array)[2])}
  \item{fp.mask.auto}{a boolean indicating if the size (in the x and y dimension) of a circular region defined around the current pixel used to search for neighbours
should grow up until the maximum amount of neighbours are found or that the x and y limits are reached}
  \item{fp.mask.step}{a numeric value indicating the step for the growth of the searching area (only if fp.mask.auto=T) }
  \item{fp.proc}{a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
"fdr" method is not implemented yet}
}

\value{returns a list containing:
'Vx' a vector containing all the indexes of the denoised pixel's neighboors
'Ix' a vector containing the denoised signal
'varx' a vector containing the variance of the denoised signal}

\author{Tiffany Lieury,Christophe Pouzat, Yves Rozenholc}






\keyword{internal}
