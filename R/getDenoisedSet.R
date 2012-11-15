getDenoisedSet <-
function
### 'getDenoisedSet' returns the denoised dataset as an array
(fp.data.array,         
### an 3D array corresponding to the dataset (the third dimension is the time)
fp.res.listdenois
### a list containing the results of the denoising procedure (cf. 'callDenoiseVoxel')
){
  fp.dim          <- dim(fp.data.array)  
  fp.data.array   <- t(apply(fp.data.array,3,c))
  fp.data.denois  <- matrix(NA,nrow(fp.data.array),ncol(fp.data.array))
  for(fp.ii in 1:length(fp.res.listdenois)){
    fp.data.denois[,fp.ii]  <- rowMeans(fp.data.array[,fp.res.listdenois[[fp.ii]]$Vx,drop=F])}
  fp.data.denois  <- array(t(fp.data.denois),dim=fp.dim)
  return(fp.data.denois)
  ###an array containing the denoised version of the dataset
}
