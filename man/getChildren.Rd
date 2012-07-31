\name{getChildren}
\alias{getChildren}
\title{getChildren}
\description{'getChildren' returns the indexes of all the direct and indirect 
neighboors of the pixel x time at index fp.idx.cur.pix
which are still inside the list of pixels x time left to clust.}
\usage{getChildren(fp.idx.cur.pix, fp.pixtoclust, fp.res.denois)}
\arguments{
  \item{fp.idx.cur.pix}{a numeric indicating the index of the current pixel x time}
  \item{fp.pixtoclust}{a vector listing the pixels x time left to clust }
  \item{fp.res.denois}{a list containing the results of the denoising procedure (cf. callDenoiseVoxel)}
}

\value{returns a vector containing the indexes of all the direct and indirect neighboors of the pixel x time at index fp.idx.cur.pix
which are still inside the list of pixels x time left to clust}

\author{Yves Rozenholc, Christophe Pouzat, Tiffany Lieury}






\keyword{internal}
