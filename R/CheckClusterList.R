CheckClusterList <-
function
### 'CheckClusterList' is used after a cluster is modified or added in the cluster list.
### It checks if there exists in the cluster list, clusters which signal is statistically homogeneous
### with the last created or modified cluster.
### If such clusters exist, the last cluster created or modified is merged with the first homogeneous cluster found.
(fp.clust.list,
### a list of clusters containing two elements:
### - Ic a matrix of dimension time x cluster containing the signal in each cluster
### - lpix a list containing the indexes of the pixels present in each cluster (as vectors) 
fp.lastchange,
### a numeric indicating the index of the cluster last modified in the object 'fp.clust.list'
fp.var,          
### a numeric indicating the variance of the dataset
fp.alpha=.05      
### a numeric value indicating the level of the statistical multitest H0
##keyword<< internal   
)
{
  #if the cluster list contains only one cluster stop and returns:
  if(length(fp.clust.list$lpix)==1) return(list(clusters=fp.clust.list,lastchange=0))
  #copies the cluster list in a temporay list
  fp.clust.listtemp       <- fp.clust.list
  #removes from clustered matrix signal the signal of the last modified cluster 
  fp.clust.listtemp$Ic    <- fp.clust.listtemp$Ic[,-fp.lastchange,drop=F]
  #removes from the pixel x time lists, the pixel x time list of the last modified cluster 
  fp.clust.listtemp$lpix  <- fp.clust.listtemp$lpix[-fp.lastchange]
  #computes the signal of the last modified cluster
  fp.Ilm     <- fp.clust.list$Ic[,fp.lastchange]
  #computes the number of members within the class
  fp.nblm    <- length(fp.clust.list$lpix[fp.lastchange])
  #computes the matrix signals in the other clusters
  fp.Ioc     <- fp.clust.listtemp$Ic
  fp.weight  <- sapply(fp.clust.listtemp$lpix,function(fp.idx) min(fp.nblm,length(fp.idx))) 
  #computes the variances 
  fp.varall  <- fp.var/(1+fp.weight)
  #returns the indexes of the clusters that are coherent with the last modified cluster
  fp.wrt     <- MultiTestH0(fp.Ioc-fp.Ilm,fp.varall,fp.alpha)
  #if there is at least one coherent cluster do 
  if(length(fp.wrt)!=0){
    #merges the last modified cluster with the first coherent cluster found
    fp.clust.listtemp$lpix[[fp.wrt[1]]] <- unique(c(fp.clust.listtemp$lpix[[fp.wrt[1]]],fp.clust.list$lpix[[fp.lastchange]]))
    #sets lastchange to a non zero value, which means this cluster will be checked until no coherence is found with other clusters
    return(list(clusters=fp.clust.listtemp,lastchange=fp.wrt[1]))}
  return(list(clusters=fp.clust.list,lastchange=0))
  ### returns a list containing:
  ### - 'clusters', the cluster list of two elements:
  ###   * Ic a matrix of dimension time x cluster containing the signal in each cluster
  ###   * lpix a list containing the indexes of the pixels present in each cluster (as vectors)
  ### - 'lastchange' the index of the cluster last modified in the cluster list 
}
