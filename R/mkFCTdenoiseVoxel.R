mkFCTdenoiseVoxel <-
function  
### 'mkFCTdenoiseVoxel' returns a function used to denoise the pixel x time at index fp.idx.pix,
### leaving out of the denoising process the pixels x time at indexes fp.idx.pix.
### This function is used to initialise the constant variables
### when the function 'callDenoiseVoxel' is called 
##keyword<< internal
(
fp.data.coord,  
### a matrix containing all the pixel x time's (x,y) coordinates
### with x in the first column and y in the second 
fp.data.matrix, 
### a matrix of pixels x time, each row corresponding to a timepoint and each column to a pixel
fp.data.matrix.norm, 
### a matrix of normalized pixels x time, each row corresponding to a timepoint and each column to a pixel 
fp.stab.var, 
### a numeric or array indicating the variance of the dataset
fp.depth=1,
### a numeric indicating the depth of a voxel
fp.alpha=.05,       
### a numeric value indicating the level of the statistical multitest H0
fp.mask.size=NULL,
### a vector indicating the size (in the x and y dimension) of a circular region defined around the current pixel used to search for neighbours
### per default fp.mask.size is defined as c(dim(fp.data.array)[1],dim(fp.data.array)[2])
fp.mask.auto=F,
### a boolean indicating if the size (in the x and y dimension) of a circular region defined around the current pixel used to search for neighbours
### should grow up until the maximum amount of neighbours are found or that the x and y limits are reached
fp.mask.step=50,
### a numeric value indicating the step for the growth of the searching area (only if fp.mask.auto=T) 
fp.proc="bonferroni"         
### a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
### "fdr" method is not implemented yet
){
 function(
   fp.idx.pix    
   ### a numeric value indicating the index of the current pixel x time to denoise
){
  ##Find the neighboors of the current pixel x time
  #sets the increasing size of the crowns with a rounded geometric serie with a0=1 and q=1.5 (+elimination of repeated values)
  fp.cercle.size  <- c(2,6,16,39,92,212,477) #increasing crown sizes used in the paper
  #maximum number of potential neighbours 
  fp.max.neighb   <- sum(fp.cercle.size)
  #we want to find the minimum size of a mask containing the most neighbours possible
  fp.goon         <- T
  fp.neighb       <- c() 
  while(length(fp.neighb)<fp.max.neighb & fp.goon == T){
    #sets the x limits of the mask, vector of TRUE/FALSE indicated whether or not x is inside the dimension 1 of the mask
    fp.xlim         <- abs(fp.data.coord[,1]-fp.data.coord[fp.idx.pix,1])<=fp.mask.size[1]
    #sets the y limits of the mask, vector of TRUE/FALSE indicated whether or not y is inside the dimension 1 of the mask
    fp.ylim         <- abs(fp.data.coord[,2]-fp.data.coord[fp.idx.pix,2])<=fp.mask.size[2]
    #gets the indexes of the coordinates inside the mask
    fp.limits                           <- fp.xlim & fp.ylim
    fp.limits[fp.neighb]   <- F        
    fp.mask                             <- which(fp.limits)
    #substract the signal at index fp.idx.pix from the raw data matrix (for pixels x time inside the mask)
    fp.diff.matrix  <- fp.data.matrix.norm[,fp.mask]-fp.data.matrix.norm[,fp.idx.pix]
    #get the indexes of fp.idx.pix's neighboors
    fp.neighb       <- c(fp.neighb,fp.mask[MultiTestH0(fp.diff.matrix,2,fp.alpha,fp.proc)])
    if(all(fp.mask.size==c(max(fp.data.coord[,1]),max(fp.data.coord[,2])))){fp.goon <- F}
    if(fp.mask.auto==F){fp.goon <- F}
    fp.mask.size[1] <- ifelse(fp.mask.size[1]+fp.mask.step<=max(fp.data.coord[,1]),fp.mask.size[1]+fp.mask.step,max(fp.data.coord[,1]))
    fp.mask.size[2] <- ifelse(fp.mask.size[2]+fp.mask.step<=max(fp.data.coord[,2]),fp.mask.size[2]+fp.mask.step,max(fp.data.coord[,2]))
    }
  if(length(fp.neighb)==1){
    return(list(Vx=fp.idx.pix,Ix=fp.data.matrix[,fp.idx.pix],varx=fp.stab.var[,fp.idx.pix]))}
  #computes the euclidean distance (in space) between fp.idx.pix and its neighboors, orders it and return their indexes in an increasing order 
  fp.dist         <- order((fp.data.coord[fp.neighb,1]-fp.data.coord[fp.idx.pix,1])^2 + (fp.data.coord[fp.neighb,2]-fp.data.coord[fp.idx.pix,2])^2+fp.depth^2*(fp.data.coord[fp.neighb,3]-fp.data.coord[fp.idx.pix,3])^2)
  #corders the neighboors in an increasing order of euclidean distance from fp.idx.pix, and adds the index of p.idx.pix first (+elimination of of repeated values)
  fp.neighb       <- unique(c(fp.idx.pix,fp.neighb[fp.dist]))
  ##Denoising procedure
  ##Vx is the neighboorhood of the pixel x time x. We want to build Vx such that is still coherent with
  ##Wx is an increasing crown around x which doesn't contain any pixel x time already in Vx
  #creates Wx, at step 1, Wx contains one pixel x time: x. fp.W is a list where Wx will be stored at each step of the algorithm
  fp.W            <- list(fp.neighb[1])
  #computes the signal Iw of the crown at step 1. fp.Iw is a matrix where Iw will be stored at each step of the algorithm (row=Iw,col=steps)
  fp.Iw           <- fp.data.matrix[,fp.neighb[1]]
  #computes the variance varw of the signal of the crown at step 1. fp.varw is a vector where varw will be stored at each step of the algorithm
  fp.varw         <- fp.stab.var[,fp.neighb[1]]
  #eliminates pixel x time x from the neighboors list
  fp.neighb       <- fp.neighb[-1]
  #creates Vx, at step 1, Vx is empty. fp.V is a list where Vx will be stored at each step of the algorithm
  fp.V            <- list(c())
  #creates the signal Iv of Vx, at step 1. fp.Iv is a matrix where Iv will be stored at each step of the algorithm (row=Iw,col=steps)
  fp.Iv           <- rep(NA,nrow(fp.data.matrix)) 
  #creates the variance varv of the signal of Vx, at step 1. fp.varv will be a matrix (time x step) where varv will be stored at each step of the algorithm
  fp.varv         <- rep(NA,nrow(fp.data.matrix)) 
  #starts at step 2 (step 1 is achieved)
  fp.myindex      <- 2
  #while there are still neighboors inside the list, and crowns to search into do
  while(length(fp.neighb)>0 & fp.myindex<=length(fp.cercle.size)){
    #returns the minimum value between the size of the crown at index fp.myindex-1 and the size of neighboors lists, ensures that one doesnt take more neighboors than available
    fp.limits           <- min(fp.cercle.size[fp.myindex-1],length(fp.neighb))
    #merges the crown's neighboors  stored at step fp.myindex-1, with the Vx's neighboors stored at step fp.myindex-1 (at step 2 Vx ={x})
    fp.V[[fp.myindex]]  <- c(fp.V[[fp.myindex-1]],fp.W[[fp.myindex-1]])
    #computes the signal Iv, average of the signals of pixels x time in Vx at step fp.myindex. Include it in the storage matrix  
    fp.Iv               <- cbind(fp.Iv,rowMeans(sapply(fp.V[[fp.myindex]],function(fp.idx) fp.data.matrix[,fp.idx])))
    #computes the variance varv, variance divided by the number of pixels x time in Vx at step fp.myindex. Include it in the storage vector 
    fp.varv             <- cbind(fp.varv,rowMeans(sapply(fp.V[[fp.myindex]],function(fp.idx) fp.stab.var[,fp.idx]))/length(fp.V[[fp.myindex]]))    
    #get the neighboors contained in the next crown
    fp.W[[fp.myindex]]  <- fp.neighb[1:fp.limits]
    #computes the signal Iw, average of the signals of pixels x time in Wx at step fp.myindex. Include it in the storage matrix  
    fp.Iw               <- cbind(fp.Iw,rowMeans(sapply(fp.W[[fp.myindex]],function(fp.idx) fp.data.matrix[,fp.idx])))
    #computes the variance varw, variance divided by the number of pixels x time in Wx at step fp.myindex. Include it in the storage vector 
    fp.varw             <- cbind(fp.varw,rowMeans(sapply(fp.W[[fp.myindex]],function(fp.idx) fp.stab.var[,fp.idx]))/length(fp.W[[fp.myindex]]))   
    #eliminates the included neighboors from the list of available neighboors
    fp.neighb           <- fp.neighb[-(1:fp.limits)]
    #tests the coherence between the signal in the crown at fp.myindex with all the previous Vx created (from step 2 since Vx={} at step 1), fp.alpha is divided by the number of Vx used at each step
    fp.testcoh          <- MultiTestH0(fp.Iv[,2:fp.myindex]-fp.Iw[,fp.myindex],fp.varv[,2:fp.myindex]+fp.varw[,fp.myindex],fp.alpha/(fp.myindex-1),fp.proc)
    #if all the Vx are not coherent with the current Wx then the denoising stops and returns Vx and Iv at index fp.myindex 
    if(length(fp.testcoh)!=length(2:fp.myindex)) return(list(Vx=unique(fp.V[[fp.myindex]]),Ix=fp.Iv[,fp.myindex],varx=fp.varv[,fp.myindex]))
    #if there are no neighboors available then the denoising stops and returns Vx and Iv at index fp.myindex 
    if(fp.limits==0)  return(list(Vx=unique(fp.V[[fp.myindex]]),Ix=fp.Iv[,fp.myindex],varx=fp.varv[,fp.myindex]))
    #next step
    fp.myindex          <- fp.myindex+1}
  #returns Vx and Iv at index fp.myindex 
  return(list(Vx=unique(fp.V[[length(fp.W)]]),Ix=fp.Iv[,length(fp.W)],varx=fp.varv[,length(fp.W)]))}
  ### returns a list containing:
  ### 'Vx' a vector containing all the indexes of the denoised pixel's neighboors
  ### 'Ix' a vector containing the denoised signal
  ### 'varx' a vector containing the variance of the denoised signal
}
