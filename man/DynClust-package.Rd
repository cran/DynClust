\name{DynClust-package}
\alias{DynClust-package}
\docType{package}
\title{non-parametric denoising and clustering method of noisy images both indexed by time and space}
\description{Two-stage method for the denoising and clustering of stack of noisy images acquired over time, based on the simple assumption that a finite sequence of noisy images both indexed by time and space is composed of noisy versions of only a limited amount of dynamic features. The aim of the method is first to denoise signals using both the spatial and temporal information contained in the data, and then cluster the denoised signals depending on their dynamic features. Two signals are considered to have similar features if their difference does not significantly deviate from zero. By comparing difference signals, no assumption is therefore made on the shape of the theoretical signals. In order for the method to be applicable to experimental data, the data should be normally distributed (or at least follow a symmetric distribution) with a constant variance. Also the number of observations n must be of the form n=d^2. Moreover, the method is based on the implicit assumption that, for a given data set, almost each dynamic feature is present in two or more pixels. The use of the method can be time consumming depending on the size of the data-array (for example about 2 hours for 130x175x60 data set on a PC  dual core processor 2.80 GHz, RAM 4 Go, Ubuntu Maverick 64 bits, R 2.12.1).}
\details{
\tabular{ll}{Package: \tab DynClust\cr
Type: \tab Package\cr
Title: \tab non-parametric denoising and clustering method of noisy images both indexed by time and space\cr
Version: \tab 1.1\cr
Date: \tab 2012-07-31\cr
Author: \tab Yves Rozenholc, Christophe Pouzat, Tiffany Lieury\cr
Maintainer: \tab Tiffany Lieury <lieuryt-dynclust@yahoo.fr>\cr
Depends: \tab R (>= 2.10)\cr
License: \tab GPL (>=2)\cr}
}
\author{Yves Rozenholc, Christophe Pouzat, Tiffany Lieury}

\keyword{ package }


