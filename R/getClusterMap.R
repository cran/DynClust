getClusterMap <-
function
### 'getClusteredMap' returns the location of each cluster as matrix
(fp.dim,         
### dimension of the dataset
fp.res.clust
### a list containing the results of the clustering procedure (cf. 'ClusteringFct')
){
  fp.dim        <- fp.dim[-length(fp.dim)]
  if(length(fp.dim)==3){
    fp.clustmap   <- array(NA,dim=fp.dim)
    fp.clustmap   <- t(apply(fp.clustmap,3,c))
  }else{
    fp.clustmap <- matrix(NA,nrow=1,ncol=prod(fp.dim))
  }    
  for(fp.ii in 1:length(fp.res.clust)){
    fp.pixcl                <- fp.res.clust[[fp.ii]]
    fp.clustmap[,fp.pixcl]  <- fp.ii}
  fp.clustmap  <- array(fp.clustmap,dim=fp.dim)
  return(fp.clustmap)
  ## returns an array containing the clustered version of the dataset
  ## returns a matrix indicating using their creation index the location of each cluster
}
