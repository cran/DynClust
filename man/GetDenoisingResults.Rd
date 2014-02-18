\name{GetDenoisingResults}
\alias{GetDenoisingResults}
\title{Get denoising step result}
\usage{
GetDenoisingResults(data.array, res.listdenois)
}
\arguments{
  \item{data.array}{a (2D or 3D)+T array containing the
  original dynamic sequence of images (the dataset). The
  last dimension is the time.}

  \item{res.listdenois}{the list resulting from the
  \code{\link{RunDenoising}} procedure applied to
  \code{data.array}. This parameter may be replaced by the
  component \code{info.den} of the former.}
}
\value{
an array with same dimension as \code{data.array}
containing the denoised version.
}
\description{
\code{GetDenoisingResults} returns the denoised version of
a dynamical image sequence as an array having the same
dimensions as the original sequence.
}
\examples{

\dontrun{
    library(DynClust)

    ## use fluorescence calcium imaging of neurons performed with Fura 2 excited at 340 nm
    data('adu340_4small',package='DynClust')

    ## Gain of the CCD camera:
    G <- 0.146
    ## readout variance of the CCD camera:
    sro2 <- (16.4)^2
    ## Stabilization of the variance to get a normalized dataset (variance=1)
    FT <- 2*sqrt(adu340_4small/G + sro2)
    FT.range = range(FT)

    ## launches the denoising step on the dataset with a statistical level of 5\%
    FT.den.tmp <- RunDenoising(FT,1,mask.size=NA,nproc=2)

    ## get the results of the denoising step
    FT.den.res <- GetDenoisingResults(FT,FT.den.tmp)

    ## plot results at time 50 in same grey scale
    par(mfrow=c(1,3))
    image(FT[,,50],zlim=FT.range,col=gray(seq(0,1,l=128)))
    title('Original')
    image(FT.den.res[,,50],zlim=FT.range,col=gray(seq(0,1,l=128)))
    title('Denoised')
    image(FT.den.res[,,50]-FT[,,50],col=gray(seq(0,1,l=128)))
    title('Residuals')
}
}
\author{
Tiffany Lieury, Christophe Pouzat, Yves Rozenholc
}
\references{
Rozenholc, Y. and Reiss, M. (2012) \emph{Preserving time
structures while denoising a dynamical image}, Mathematical
Methods for Signal and Image Analysis and Representation
(Chapter 12), Florack, L. and Duits, R. and Jongbloed, G.
and van~Lieshout, M.-C. and Davies, L. Ed.,
Springer-Verlag, Berlin

Lieury, T. and Pouzat, C. and Rozenholc, Y. (submitted)
\emph{Spatial denoising and clustering of dynamical image
sequence: application to DCE imaging in medicine and
calcium imaging in neurons}
}
\seealso{
\code{\link{RunDenoising}}
}

