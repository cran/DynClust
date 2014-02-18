\name{GetClusteringResults}
\alias{GetClusteringResults}
\title{Get clustering step result}
\usage{
GetClusteringResults(data.array, res.listdenois, res.cluster)
}
\arguments{
  \item{data.array}{a (2D or 3D)+T array containing the
  original dynamic sequence of images (the dataset). The
  last dimension is the time.}

  \item{res.listdenois}{the list resulting from the
  \code{\link{RunDenoising}} procedure applied to
  \code{data.array}. This parameter may be replaced by the
  component \code{info.den} of the former.}

  \item{res.cluster}{the list resulting from a call to
  \code{\link{RunClustering}}}
}
\value{
a list containing two components \code{clust.array} and
\code{clust.map}. \code{clust.array} is an array with same
dimension as the original sequence \code{data.array}
containing the clustered version. \code{clust.map} is an
array with only spatial dimensions of \code{data.array}
whose elements provide the cluster number at each location.
}
\description{
\code{GetClusteringResults} returns the results of the
clustering procedure \code{RunClustering}
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
    
    ## launches the clustering step on the dataset with a statistical level of 5\%
    FT.clust.tmp  <- RunClustering(FT,FT.den.tmp,nproc=2)
    n.cluster <- length(FT.clust.tmp$clusters)
    print(paste(n.cluster,'clusters using variance set to',sqrt(FT.den.tmp$var),'^2'))
    
    ## get the classified version of the data array and the map of the clusters
    FT.clust.res <- GetClusteringResults(FT,FT.den.tmp,FT.clust.tmp)
    
    ## plotting results of the clusterization
    par(mfrow=c(2,2))
    image(FT.clust.res$clust.map,col=rainbow(n.cluster))
    title('Cluster map')
    matplot(FT.clust.res$clust.center,col=rainbow(n.cluster),type="l",lwd=0.1,lty=1)
    title('Cluster centers')

    ## and more: original and clustered slices at time 50
    image(FT[,,50],zlim=FT.range,col=grey(seq(0,1,length=n.cluster)))
    title('Original sequence at time 50')
    image(FT.clust.res$clust.array[,,50],zlim=FT.range,col=grey(seq(0,1,length=n.cluster)))
    title('Clustered sequence at time 50')

    ####################################################################################
    ## reapply clustering with twice the nominal variance: forces stronger clustering ##
    ####################################################################################

    ## launches the denoising step on the dataset with a statistical level of 5\%
    FT.den.tmp <- RunDenoising(FT,2,mask.size=NA,nproc=2)

    ## launches the clustering step on the dataset with a statistical level of 5\%
    FT.clust.tmp  <- RunClustering(FT,FT.den.tmp,nproc=2)
    n.cluster <- length(FT.clust.tmp$clusters)
    print(paste(n.cluster,'clusters using twice the nominal variance'))
    
    ## get the classified version of the data array and the map of the clusters
    FT.clust.res <- GetClusteringResults(FT,FT.den.tmp,FT.clust.tmp)
    
    ## plotting results of the clusterization
    par(mfrow=c(2,2))
    image(FT.clust.res$clust.map,col=rainbow(n.cluster))
    title('Cluster map')
    matplot(FT.clust.res$clust.center,col=rainbow(n.cluster),type="l",lwd=0.1,lty=1)
    title('Cluster centers')

    ## and more: original and clustered slices at time 50
    image(FT[,,50],zlim=FT.range,col=grey(seq(0,1,length=n.cluster)))
    title('Original sequence at time 50')
    image(FT.clust.res$clust.array[,,50],zlim=FT.range,col=grey(seq(0,1,length=n.cluster)))
    title('Clustered sequence at time 50')
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
\code{\link{RunDenoising}}, \code{\link{RunClustering}}
}

