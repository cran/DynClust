\name{mkFCTdenoiseVoxel}
\alias{mkFCTdenoiseVoxel}
\title{mkFCTdenoiseVoxel}
\description{'mkFCTdenoiseVoxel' returns a function used to denoise the pixel x time at index fp.idx.pix,
leaving out of the denoising process the pixels x time at indexes fp.idx.pix.
This function is used to initialise the constant variables
when the function 'callDenoiseVoxel' is called }
\usage{mkFCTdenoiseVoxel(fp.data.coord, fp.data.matrix, fp.stab.var, 
    fp.alpha = 0.05, fp.proc = "bonferroni")}
\arguments{
  \item{fp.data.coord}{a matrix containing all the pixel x time's (x,y) coordinates
with x in the first column and y in the second }
  \item{fp.data.matrix}{a matrix of pixels x time, each row corresponding to a timepoint and each column to a pixel}
  \item{fp.stab.var}{a numeric or array indicating the variance of the dataset}
  \item{fp.alpha}{a numeric value indicating the level of the statistical multitest H0}
  \item{fp.proc}{a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
"fdr" method is not implemented yet}
}

\value{returns a list containing:
'Vx' a vector containing all the indexes of the denoised pixel's neighboors
'Ix' a vector containing the denoised signal
'varx' a vector containing the variance of the denoised signal}

\author{Yves Rozenholc, Christophe Pouzat, Tiffany Lieury}






\keyword{internal}
