\name{RunDenoising}
\alias{RunDenoising}
\title{Denoising step of a dynamical image sequence}
\usage{
RunDenoising(data.array, data.var = 1, depth = 1, alpha = 0.05,
  mask.size = NA, nproc = 1, enhStart = ifelse(is.null(var), 2, 1))
}
\arguments{
  \item{data.array}{a (2D or 3D)+T array containing the
  dynamic sequence of images (the dataset). The last
  dimension is the time.}

  \item{data.var}{a numeric indicating the variance of the
  dataset (default 1). If set to NULL, the variance is
  computed using a baseline image. See \code{enhStart}
  parameter.}

  \item{depth}{a numeric indicating the depth of a voxel.}

  \item{alpha}{a numeric value indicating the global level
  of the multitest.}

  \item{mask.size}{a vector indicating the size of the
  spatial hypercube defined around voxels used to search
  for neighbors.

  If NA (default):
  sqrt(dim(data.array)[1:length(dim(data.array))-1]).

  If NULL (complete image):
  dim(data.array)[1:length(dim(data.array))-1]}

  \item{nproc}{a numeric value indicating the number of
  processors used for parallel computation.}

  \item{enhStart}{an integer, if larger than 1, a baseline
  is computed as a median image obtain from time indexes
  between 1 and enhStart-1. Default value
  \code{ifelse(is.null(var),2,1)}.}
}
\value{
a list containing: \itemize{ \item \code{info.den}, a list
of list whose length is the number of voxels, each sub-list
contains the result of buildEstimate for one voxel. \item
\code{data.proj}, the projections of the dynamics. a list
containing a denoised version of the dataset as an array,
as well as a list for which each element contains a list
with the voxel index, the indexes of its neighboors, the
resulting denoised signal, and the variance of the denoised
signal \item \code{var}, a numeric providing the known
variance }
}
\description{
Performs the denoising step of a dynamic sequence of
images. It is also the first step of the clustering.

The denoising procedure is iteratively applied on each
voxel.

The denoised version of Fx is obtained with a three stages
procedure: 1) Selection of time-homogeneous voxels in the
sub-mask around the voxel of interest; 2) Growth of spatial
neighborhoods made of time-homogeneous voxels obtained at
stage 1 with sizes growing geometrically---each neiborhood
is associated to a denoised version by averaging over its
members; 3) Selection of the largest spatial neighborhood
such that its associated denoised version is
time-homogeneous with all the previous ones.

Time homogeneity is tested with function
\code{MultiTestH0}.

Further details about the denoising method and the
statistical test of homogeneity can be found in the
references.
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
\code{\link{GetDenoisingResults}},
\code{\link{MultiTestH0}}.
}

