getClusteredSet <-
function
### 'getClusteredSet' returns the dataset where each pixel x time was replaced by its center
(fp.data.array,         
### an 3D array corresponding to the dataset (the third dimension is the time)
fp.res.denois,
### a list containing the results of the denoising procedure (cf. 'callDenoiseVoxel')
fp.res.clust
### a list containing the results of the clustering procedure (cf. 'ClusteringFct')
){
  fp.dim        <- dim(fp.data.array)
  fp.data.array <- getDenoisedSet(fp.data.array,fp.res.denois)
  fp.data.array <- t(apply(fp.data.array,3,c))
  fp.data.clust <- matrix(NA,nrow(fp.data.array),ncol(fp.data.array))
  for(fp.ii in 1:length(fp.res.clust)){
    fp.pixcl  <- fp.res.clust[[fp.ii]]
    fp.data.clust[,fp.pixcl]  <- rowMeans(fp.data.array[,fp.pixcl,drop=F])}
  fp.data.clust <- array(t(fp.data.clust),dim=fp.dim)
  return(fp.data.clust)
  ## returns an array containing the clustered version of the dataset
}
