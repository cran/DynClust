callDenoiseVoxel <-
structure(function
### 'callDenoiseVoxel' performs the denoising step of the clustering method.
### The function returns a list containing a denoised version of the dataset as an array,
### as well as a list for which each element contains a list with the pixel index,
### the indexes of its neighboors, the resulting denoised signal, and the variance of the denoised signal   
(fp.data.array,         
### an 3D array corresponding to the dataset (the third dimension is the time) 
fp.stab.var,           
### a numeric or array indicating the variance of the dataset
fp.alpha=.05,       
### a numeric value indicating the level of the statistical multitest H0
fp.proc="bonferroni"         
### a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
### "fdr" method is not implemented yet 
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
