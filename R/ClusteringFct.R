ClusteringFct <-
structure(function
### 'ClusteringFct' clusters together the denoised pixels x time
### with  statistically homogeneous signals at a level 'fp.alpha'
(fp.res.listdenois,
### a list containing the results of the denoising procedure (cf. 'callDenoiseVoxel')
fp.alpha=.05,       
### a numeric value indicating the level of the statistical multitest H0
fp.proc="bonferroni"         
### a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
### "fdr" method is not implemented yet
##details<< As for the denoising procedure, the classes are built such that within each class
##details<< the clustered signals are statistically coherent. The clustering procedure is a recursive
##details<< algorithm, and the main loop of which can be decomposed into the following steps:
##details<< - the list of all the pixels, Lv, is ordered by increasing neighborhood sizes;
##details<< - a class is defined from a single pixel: at step 1, the class is defined using
##details<< the first pixel of Lv. The coherence between the signal in the current pixel and
##details<< all other pixels in Lv is statistically tested. The pixels for which the null hypothesis
##details<< was not rejected are removed from Lv and included in the current cluster. The estimated
##details<< signal of the cluster is computed as the averaged signal over all the pixels of the cluster.
##details<< - the estimated signal of the new class is statistically tested for coherence with the
##details<< estimated signal of all other existing clusters;
##details<< - statistically coherent clusters are merged.
##details<< The algorithm stops when Lv is empty and when the clusters cannot be merged anymore.
##details<< The number of clusters could be decreased by tuning the level alpha or the variance
##details<< during the denoising procedure.
##details<< However, by tuning one of these parameters, users expose themselves to the risk of 
##details<< obtaining under or over-smoothed estimated signals in the clusters.
##details<< Details about the statistical test for coherence can be found in the references.
##references<< Rozenholc, Y. and Reiss, M. (2012) _Preserving time structures while denoising a dynamical image_,
##references<< Mathematical Methods for Signal and Image Analysis and Representation (Chapter 12), 
##references<< Florack, L. and Duits, R. and Jongbloed, G. and van~Lieshout, M.-C. and Davies, L.
##references<< Ed., Springer-Verlag, Berlin
){
  #computes the sizes of every pixels x time neighboorhood
  fp.res.denois         <- fp.res.listdenois$details
  fp.dimarray           <- dim(fp.res.listdenois$denois3D)
  fp.neighlen           <- sapply(fp.res.denois,function(fp.idx) length(fp.res.denois$Vx))
  #orders the sizes in a decreasing order will become the list of pixel x time to clust
  fp.neighlen.sort.idx  <- sort(fp.neighlen,decreasing=T,index.return=T)$ix
  #creates a function for clustering, initialises the constant variables
  fp.ClusteringInnerFct <- mkClusteringInnerFct(fp.res.denois,fp.alpha,fp.proc)
  #step 1 of the clustering, where the list of clusters = NULL
  fp.res.clust          <- fp.ClusteringInnerFct(fp.neighlen.sort.idx[1],0,NULL,fp.neighlen.sort.idx)
  #while the list of pixels x time to clust is not empty do
  fp.ncolmat          <- length(fp.res.clust$pixtoclust)
  fp.prog             <- follow_progress(fp.idx.curr=fp.ncolmat,fp.max.iter=fp.ncolmat,fp.prog=seq(0,100,by=1))
  while(length(fp.res.clust$pixtoclust)>0){
    #continues the clustering updating at each step the list of pixel x time to clust, and choosing the next pixel x time to clust
    fp.res.clust  <- fp.ClusteringInnerFct(fp.res.clust$pixtoclust[1],fp.res.clust$lastchange,fp.res.clust$clusters,fp.res.clust$pixtoclust)
    fp.prog       <- follow_progress(length(fp.res.clust$pixtoclust),fp.ncolmat,fp.prog)
  }
  
  #sorts the cluster by their neighboorhood size in the decreasing order
  fp.order                      <- order(sapply(fp.res.clust$clusters$lpix,length),decreasing=T)
  fp.res.clust$clusters$lpix    <- lapply(fp.order,function(fp.idx){fp.res.clust$clusters$lpix[[fp.idx]]}) 
  fp.res.clust$clusters$Ic      <- fp.res.clust$clusters$Ic[,fp.order,drop=F]
  fp.res.clust$clusters$varc    <- fp.res.clust$clusters$varc[,fp.order,drop=F]
  
  #constructs a 3D data array with the clustered signals
  fp.clust.mat       <- matrix(NA,nrow=fp.dimarray[3],ncol=fp.dimarray[1]*fp.dimarray[2])
  for(ii in 1:length(fp.res.clust$cluster$lpix))
    fp.clust.mat[,fp.res.clust$clusters$lpix[[ii]]]  <- fp.res.clust$clusters$Ic[,ii]
  fp.clust.array    <- array(t(fp.clust.mat),dim=fp.dimarray)
  
    
  #constructs a 3D data array with the variance of the clustered signals
  fp.clust.mat       <- matrix(NA,nrow=fp.dimarray[3],ncol=fp.dimarray[1]*fp.dimarray[2])
  for(ii in 1:length(fp.res.clust$cluster$lpix))
    fp.clust.mat[,fp.res.clust$clusters$lpix[[ii]]]  <- fp.res.clust$clusters$varc[,ii]
  fp.clust.array.var  <- array(t(fp.clust.mat),dim=fp.dimarray)

  #constructs a map of the location of the clusters
  fp.img.clust <- rep(NA,fp.dimarray[1]*fp.dimarray[2])
  for(ii in 1:length(fp.res.clust$cluster$lpix)){
    fp.img.clust[fp.res.clust$clusters$lpix[[ii]]]  <- ii}
  fp.img.clust <- matrix(fp.img.clust,fp.dimarray[1],fp.dimarray[2])
  
  return(list(details=fp.res.clust$clusters,clust3D=fp.clust.array,clust3Dvar=fp.clust.array.var,clustmap=fp.img.clust))
  ### returns a list containing:
  ### - 'details' a list containing for each cluster a list of 3 objects:
  ###   * Ic a matrix of dimension time x cluster containing the signal in each cluster
  ###   * varc a matrix of dimension time x cluster containing variance of the signal in each cluster
  ###   * lpix a list containing the indexes of the pixels present in each cluster (as vectors)
  ### - 'clust3D' an array containing the clustered version of the dataset
  ### - 'clustmap' a matrix indicating using their creation index the location of each cluster
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
# #launches the clustering step on the dataset with a statistical level of 5%
# clustres      <- ClusteringFct(denoisres,.05)
# x11()
# matplot(clustres$details$Ic,t="l",lty=1,lwd=2,bty="n")
# x11()
# par(mar=rep(0,4))
# image(clustres$clustmap,col=rainbow(max(clustres$clustmap)))
# x11()
# par(mar=rep(0,4))
# image(apply(adu340_4small,1:2,mean),col=grey(seq(0,1,len=2^8)))
})
