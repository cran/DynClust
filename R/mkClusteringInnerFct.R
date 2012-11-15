mkClusteringInnerFct <-
function
### 'mkClusteringInnerFct' returns the clustering function
(fp.res.denois,
### resulting object obtained after denoising
fp.var,          
### a numeric indicating the variance of the dataset
fp.alpha=.05      
### a numeric value indicating the level of the statistical multitest H0
##keyword<< internal  
){
  function(
    fp.idx.cur.pix,
    ### index of the current pixel x time x
    fp.lastchange,
    ### index of the cluster 
    fp.clust.list,
    ### object containing the result of clustering
    fp.pixtoclust
    ### list of the pixel x time left to clust
    ){
      #at step 1, the cluster list is null, initialises it, list object with Ic matrix of signals Ic inside each cluster, and list of the Neighboors for each cluster in lpix
      if(is.null(fp.clust.list)) fp.clust.list  <- list(Ic=c(),lpix=list())
      #if the last change occurred at index 0 do
      if(fp.lastchange==0){
        #gets the complete neighboors list to evaluate for clusting (see getChildren)
        fp.newcluster <- fp.pixtoclust       
        #computes the matrix of denoised signals of every neighboor in fp.newcluster
        fp.Idm        <- sapply(fp.newcluster,function(idx) fp.res.denois[[idx]]$Ix)
        #computes the denoised signal of x
        fp.Idx        <- fp.res.denois[[fp.idx.cur.pix]]$Ix
        #computes the variances for each denoised signal
        fp.vardm      <- fp.var/sapply(fp.newcluster,function(fp.idx) length(fp.res.denois[[fp.idx]]$Vx))
        #computes the variance for the denoised signal at x
        fp.vardx      <- fp.var/length(fp.res.denois[[fp.idx.cur.pix]]$Vx)
        fp.newvar     <- fp.vardm+fp.vardx   
        #tests the coherence between the denoised signal at x, and the denoised signals proposed as potential neighboors 
        #eliminates the pixels x time that are not coherent with x
        fp.newcluster <- fp.newcluster[MultiTestH0(fp.Idm-fp.Idx,fp.newvar,fp.alpha)]
        #one cluster has been created, the size of the cluster list is incremented
        fp.lastchange <- length(fp.clust.list$lpix)+1
      }else{
        #sets the neighboors list to the list contained in the cluster's pixels x time list at index fp.lastchange
        fp.newcluster <- fp.clust.list$lpix[[fp.lastchange]]}
      #computes a signal to represent the average element of the class
      fp.newclustcenter <- list(lpix=fp.newcluster,
                                Ic  = rowMeans(sapply(fp.newcluster,function(idx) fp.res.denois[[idx]]$Ix)))
      #updates the list of pixels x time left to clust
      fp.pixtoclust     <- updateList(fp.pixtoclust,fp.newcluster)
      #if no new cluster were created do (or when the first class is created)
      if(fp.lastchange<=length(fp.clust.list$lpix)){
        #updates the cluster's signal
        fp.clust.list$Ic[,fp.lastchange]   <- fp.newclustcenter$Ic}
      else{
        #if a new cluster is created, had the new cluster's signal to matrix of clusters' signals
        fp.clust.list$Ic   <- cbind(fp.clust.list$Ic,fp.newclustcenter$Ic)}
      #updates the list of pixels x time in the visited cluster whether it is a new cluster or not  
      fp.clust.list$lpix[[fp.lastchange]]     <- fp.newclustcenter$lpix
      #checks if the visiter cluster is or not coherent with other existing cluster, returns the cluster's index of the last change 
      fp.checkcl     <- CheckClusterList(fp.clust.list,fp.lastchange,fp.var,fp.alpha)
      #update the index of the last change 
      fp.lastchange  <- fp.checkcl$lastchange
      #update the clusters list
      fp.clust.list  <- fp.checkcl$cluster
      #returns the results of clustering
      return(list(clusters=fp.clust.list,lastchange=fp.lastchange,pixtoclust=fp.pixtoclust))
      ### returns a list containing:
      ### - 'clusters', the cluster list of two elements:
      ###   * Ic a matrix of dimension time x cluster containing the signal in each cluster cluster
      ###   * lpix a list containing the indexes of the pixels present in each cluster (as vectors)
      ### - 'lastchange' the index of the cluster last modified in the cluster list
      ### - 'pixtoclust' a vector containing the indexes of the pixels left for clustering    
    }
}
