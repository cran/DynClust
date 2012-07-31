callDenoiseVoxel <-
structure(function
### 'callDenoiseVoxel' performs the denoising step of the clustering method.
### The function returns a list containing a denoised version of the dataset as an array,
### as well as a list for which each element contains a list with the pixel index,
### the indexes of its neighboors, the resulting denoised signal, and the variance of the denoised signal   
(fp.data.array,         
### an 3D array corresponding to the dataset (the third dimension is the time).
### The number of observations n must be of the form n=2^d
fp.stab.var,           
### a numeric or array indicating the variance of the dataset
fp.alpha=.05,       
### a numeric value indicating the level of the statistical multitest H0
fp.proc="bonferroni"         
### a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
### "fdr" method is not implemented yet 
##details<< The denoising procedure is applied on each signal. Given a pixel at location x, 
##details<< the difference 3D data-array is constructed by subtracting the signal at location x, Fx,
##details<< from all the other signals of the data set. The time-homogeneity is statistically tested
##details<< for each difference signal at a level alpha. The locations of the pixels for which the null
##details<< hypothesis was not rejected are saved, and will be listed as neighbors of the pixel x.
##details<< These neighbors are then used to construct a denoised estimation of Fx by averaging the
##details<< closest (in space) time-homogeneous signals selected among all the listed neighbors.
##details<< Since it is expected that the closer a neighbor is from the current pixel, the most likely
##details<< their signals will be coherent, neighbors are ordered by proximity to the current pixel.
##details<< Several estimations of the denoised Fx are constructed by averaging the signals of an
##details<< increasing number of neighbors. The size of the neighborhood is increasing geometrically
##details<< to optimize the computational costs while still obtaining estimates achieving a good statistical
##details<< trade-off between bias and variance in step.
##details<< The estimation of the denoised Fx using the largest neighborhood size (hence the highest denoising level),
##details<< is then compared to all other estimations. If all the estimations are statistically coherent with each other,
##details<< the estimation of the denoised signal located at pixel x using the largest neighborhood size becomes
##details<< the final denoised Fx. However, if at least one of the estimations is not statistically coherent with
##details<< the one having the largest neighborhood size, the estimation of the denoised signal located at pixel x
##details<< using the largest neighborhood size is eliminated from the possible estimations.
##details<< The process is repeated until all estimates are statistically close to the estimation
##details<< computed with the largest neighborhood size.
##details<< The number of clusters could be decreased by tuning the level alpha or the variance
##details<< during the denoising procedure.
##details<< However, by tuning one of these parameters, users expose themselves to the risk of 
##details<< obtaining under or over-smoothed estimated signals in the clusters. 
##details<< Details about the denoising method and the statistical test for coherence can be found in the references.
##references<< Rozenholc, Y. and Reiss, M. (2012) _Preserving time structures while denoising a dynamical image_,
##references<< Mathematical Methods for Signal and Image Analysis and Representation (Chapter 12), 
##references<< Florack, L. and Duits, R. and Jongbloed, G. and van~Lieshout, M.-C. and Davies, L.
##references<< Ed., Springer-Verlag, Berlin
){
  fp.n  <- dim(fp.data.array)[3] 
  if(log2(fp.n) != round(log2(fp.n))){stop("the time dimension is not of the form n=2^d => log2(n) is not an integer")}
  if(!is.numeric(fp.stab.var) & !is.array(fp.stab.var)){stop("fp.stab.var should be an array or a numeric")}
  if(is.array(fp.stab.var) & !all(dim(fp.stab.var)== dim(fp.data.array))){stop("dim(fp.stab.var) and dim(fp.data.array) are different")}
  if(is.numeric(fp.stab.var)) fp.stab.var  <- array(fp.stab.var,dim=dim(fp.data.array)) 

  #transforms the 3 array into a 2D matrix with nrow = time, ncol pixels x time
  fp.data.matrix      <- t(apply(fp.data.array,3,c))
  #transforms fp.stab.var into the appropriate object (into a 2D matrix with nrow = time, ncol pixels x time)
  fp.stab.var         <- t(apply(fp.stab.var,3,c))
  #matrix of all the 3D coordinates at column indexes of fp.data.matrix corresponds to row indexes in fp.data.coord
  fp.data.coord       <- cbind(x=rep(1:nrow(fp.data.array),ncol(fp.data.array)),y=rep(1:ncol(fp.data.array),each=nrow(fp.data.array)))
  #creates the function fp.fctDenoiseVoxel and initialises parameters that will not change until the end of the analysis
  fp.DenoiseVoxel     <- mkFCTdenoiseVoxel(fp.data.coord,fp.data.matrix,fp.stab.var,fp.alpha,fp.proc)
  #creates a list of size number of pixels x time, where the results will be stored
  fp.res.list         <- vector("list",ncol(fp.data.matrix))
  #for each pixel x time do
  fp.ncolmat          <- ncol(fp.data.matrix)
  fp.prog             <- follow_progress()
  for(fp.idx.pix in 1:fp.ncolmat){
    #creates the function fp.fctDenoise and initialises the parameters 
    fp.res.list[[fp.idx.pix]] <- fp.DenoiseVoxel(fp.idx.pix)
    fp.prog                   <- follow_progress(fp.idx.curr=fp.idx.pix,fp.max.iter=fp.ncolmat,fp.prog=fp.prog)
  }
  fp.denois.array  <- array(t(sapply(1:length(fp.res.list),function(idx){fp.res.list[[idx]]$Ix})),dim=dim(fp.data.array))
  list(details=fp.res.list,denois3D=fp.denois.array)
  ### returns a list containing:
  ### 'details' a list containing for each pixel x time a list of 4 objects:
  ### - 'Vx' a vector containing all the indexes of the denoised pixel's neighboors
  ### - 'Ix' a vector containing the denoised signal
  ### - 'varx' a vector containing the variance of the denoised signal
  ### 'denois3D' an array containing the denoised version of the dataset
}, ex = function(){
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
# #launches the denoising step on the dataset with a statistical level of 5%
# denoisres     <- callDenoiseVoxel(adu340_4small,FT_varhat,.05,"bonferroni")
# #computes the average over time for each pixelxtime before denoising
# avg_before    <- apply(adu340_4small,1:2,mean)
# #computes the average over time for each pixelxtime after denoising
# avg_after     <- apply(denoisres$denois3D,1:2,mean)
# #plotting results
# image(1:nrow(avg_before),1:ncol(avg_before),avg_before)
# x11();image(1:nrow(avg_after),1:ncol(avg_after),avg_after)
})
