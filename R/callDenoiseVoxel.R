callDenoiseVoxel <-
structure(function
### 'callDenoiseVoxel' performs the denoising step of the clustering method.
### The function returns a list containing a denoised version of the dataset as an array,
### as well as a list for which each element contains a list with the pixel index,
### the indexes of its neighboors, the resulting denoised signal, and the variance of the denoised signal   
(fp.data.array,         
### an 3D array corresponding to the dataset (the third dimension is the time).
### The number of observations n must be of the form n=2^d
fp.var=1,          
### a numeric indicating the variance of the dataset
fp.depth=1,
### a numeric indicating the depth of a voxel
fp.alpha=.05,       
### a numeric value indicating the level of the statistical multitest H0
fp.mask.size=sqrt(dim(fp.data.array)[1:2]),
### a vector indicating the size (in the x and y dimension) of a rectangular region defined around the current pixel used to search for neighbours
### per default fp.mask.size = sqrt(dim(fp.data.arra)[1:2]), set it to fp.mask.size = NULL to use the whole dataset
fp.nproc=1
### a numeric value indicating the number of processors to run parallel computations
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
  ############################# DIM INITIALIZATION #############################
  #data dimensions
  fp.dim              <- dim(fp.data.array)
  #dimension variables
  fp.nx               <- fp.dim[1]
  fp.ny               <- fp.dim[2]
  fp.n                <- fp.dim[length(fp.dim)]
  fp.nz               <- ifelse(length(fp.dim)==4,fp.dim[3],1)
  ############################# TESTING FEASIBILITY #############################
  #tests can the data be analyzed
  if(fp.nproc<1){stop("fp.nproc must be at least equal to 1")}
  if(log2(fp.n) != round(log2(fp.n))){stop("the time dimension is not of the form n=2^d => log2(n) is not an integer")}
  if(!is.null(dim(fp.var)) | length(fp.var)>1){stop("fp.var should be a numeric")}
  if(length(fp.dim)!=4 & length(fp.dim)!=3){stop("number of dimension must be either 3 or 4")}
  if(length(fp.dim)==3){ dim(fp.data.array) <- fp.dim <- c(fp.nx,fp.ny,fp.nz,fp.n)}
  if(!is.null(fp.mask.size))  fp.mask.size  <- round(abs(fp.mask.size))
  ############################# CREATING LOCAL FUNCTIONS #############################
  # creates the functions needed for denoising and initialises static parameters
  ## local function used to compute for the projections of all the dyadic partitions for each voxel
  ## the function returns a matrix of ncol voxels and nrow projections
  RetProj <- function(){
    fp.nrow <- fp.n
    fp.ncol <- (fp.nx*fp.ny*fp.nz)
    fp.iter <- log2(fp.nrow)-1
    fp.to   <- 0
    #projection matrix
    fp.proj <- matrix(NA,nrow=fp.nrow-2,ncol=fp.ncol)
    for(fp.idx in 1:fp.iter){
      fp.from                 <- fp.to+1
      fp.to                   <- fp.from+fp.nrow/2-1     
      fp.proj[fp.from:fp.to,] <- (fp.data.array[seq(1,fp.nrow,by=2),,drop=F]+fp.data.array[seq(2,fp.nrow,by=2),,drop=F])/2
      fp.data.array          <- fp.proj[fp.from:fp.to,,drop=F]
      fp.nrow                 <- nrow(fp.data.array)
    }
    return(fp.proj)
  }  

  ## local function used to find the neighboors of the current pixel (at index fp.idx.pix),
  ## here neighboors are defined as the pixels that are time-homogeneous with the current one     
  findNeighboors <- function(fp.idx.pix){
    ##Find the neighboors of the current pixel x time
    if(!is.null(fp.mask.size)){    
      fp.xliminf    <- max(c(1,fp.data.coord[fp.idx.pix,1]-fp.mask.size[1]/2))
      fp.xlimsup    <- min(c(fp.nx,fp.data.coord[fp.idx.pix,1]+fp.mask.size[1]/2))
      fp.yliminf    <- max(c(1,fp.data.coord[fp.idx.pix,2]-fp.mask.size[2]/2))
      fp.ylimsup    <- min(c(fp.ny,fp.data.coord[fp.idx.pix,2]+fp.mask.size[2]/2))
      fp.xyz        <- fp.arr.coord[fp.xliminf:fp.xlimsup,fp.yliminf:fp.ylimsup,]
      fp.mask       <- which(fp.data.coord[,"idx"] %in% fp.xyz)
      #substract the signal at index fp.idx.pix from the raw data matrix (for pixels x time inside the mask)      
      #get the indexes of fp.idx.pix's neighboors
      fp.neighb       <- fp.mask[MultiTestH0(fp.data.array[,fp.mask]-fp.data.array[,fp.idx.pix],2*fp.var,fp.alpha)]
      return(fp.neighb)
    }  
    #substract the signal at index fp.idx.pix from the raw data matrix
    #get the indexes of fp.idx.pix's neighboors
    fp.neighb       <- MultiTestH0(fp.data.array-fp.data.array[,fp.idx.pix],2*fp.var,fp.alpha)
    return(fp.neighb)
    }  

  ## local function used to return the closest time-homegeneous neighboors to the current pixel
  ## take as parameters a numeric value indicating the index of the current pixel x time to denoise
  ## and a vector containing all the indexes of the denoised pixel's neighboors
  returnClosests <- function(fp.idx.pix,fp.neighb){    
      fp.coord.pix    <- fp.data.coord[fp.idx.pix,]
      #computes the euclidean distance (in space) between fp.idx.pix and its neighboors, orders it and return their indexes in an increasing order 
      fp.dx           <- (fp.data.coord[fp.neighb,1]-fp.coord.pix[1])
      fp.dy           <- (fp.data.coord[fp.neighb,2]-fp.coord.pix[2])
      fp.dz           <- fp.depth*(fp.data.coord[fp.neighb,3]-fp.coord.pix[3])
      fp.dist         <- order(fp.dx^2+fp.dy^2+fp.dz^2)
      return(unique(c(fp.idx.pix,fp.neighb[fp.dist])))
      ### returns a vector containing all the indexes of the denoised pixel's neighboors
      }  

  ## local function returning a the denoised estimation of the time dynamics and the current pixel neighboors
  ## used to compute it, the function takes as parameters
  ## a numeric value indicating the index of the current pixel x time to denoise
  ## a vector containing all the indexes of the denoised pixel's neighboors
  ## a matrix of pixels x time from each neighboors found, each row corresponding to a timepoint and each column to a pixel
  buildEstimate <-  function(fp.idx.pix, fp.neighb){
    fp.neighboors.matrix  <- fp.data.array[,fp.neighb]
    if(length(fp.neighb)==1)
      return(list(Vx=fp.idx.pix,Ix=fp.neighboors.matrix))
    ##Denoising procedure
    ##V1 à V8 are the neighboorhoods of the pixel x time x of increasing size. Here we want to find the biggest Vx 
    ##which is still coherent with the smaller ones.
    #sets the increasing size of the crowns with a rounded geometric serie with a0=1 and q=1.5 (+elimination of repeated values)
    fp.cercle.size  <- c(1,2,6,16,39,92,212,477) #increasing crown sizes used in the paper
    fp.cercle.l     <- length(fp.cercle.size)
    #creates Vx, at step 1, Vx is empty. fp.V is a list where Vx will be stored at each step of the algorithm
    fp.V            <- list(c())
    #creates the signal Iv of Vx, at step 1. fp.Iv is a matrix where Iv will be stored at each step of the algorithm (row=Iw,col=steps)
    fp.Iv           <- matrix(NA,nrow(fp.neighboors.matrix),ncol=fp.cercle.l)
    #creates the variance varv of the signal of Vx, at step 1. fp.varv will be a matrix (time x step) where varv will be stored at each step of the algorithm
    fp.varv         <- rep(NA,ncol=fp.cercle.l)
    for(fp.myindex in 1:fp.cercle.l){
      fp.limits           <- min(sum(fp.cercle.size[1:fp.myindex]),length(fp.neighb))
      fp.V[[fp.myindex]]  <- fp.neighb[1:fp.limits]
      fp.Iv[,fp.myindex]  <- rowMeans(fp.neighboors.matrix[,1:fp.limits,drop=F])
      fp.varv[fp.myindex] <- fp.var/length(fp.V[[fp.myindex]])
    }
   
    while(fp.myindex>=2){    
      fp.testcoh          <- MultiTestH0(fp.Iv[,1:(fp.myindex-1)]-fp.Iv[,fp.myindex],fp.varv[1:(fp.myindex-1)]+fp.varv[fp.myindex],fp.alpha/(fp.myindex-1))  
      if(length(fp.testcoh)==(fp.myindex-1)){
        return(list(Vx=unique(fp.V[[fp.myindex]]),Ix=fp.Iv[,fp.myindex]))
      }else{
        fp.myindex  <- fp.myindex-1
      }}
      return(list(Vx=unique(fp.V[[fp.myindex]]),Ix=fp.Iv[,fp.myindex])) 
    ### returns a list containing:
    ### 'Vx' a vector containing all the indexes of the denoised pixel's neighboors
    ### 'Ix' a vector containing the denoised signal
    }

  ############################# DATA REFORMATING #############################  
  #matrix of all the 3D coordinates at column indexes of fp.data.matrix corresponds to row indexes in fp.data.coord
  fp.data.coord       <- cbind(x=rep(1:fp.nx,len=fp.nx*fp.ny*fp.nz),y=rep(rep(1:fp.ny,each=fp.nx),len=fp.nx*fp.ny*fp.nz),z=rep(rep(1:fp.nz,each=fp.ny*fp.nx),len=fp.nx*fp.ny*fp.nz),idx=1:(fp.nx*fp.ny*fp.nz))
  fp.arr.coord        <- array(1:(fp.nx*fp.ny*fp.nz),dim=c(fp.nx,fp.ny,fp.nz))
  #transforms the 3D or 4D data set into a matrix
  fp.data.array       <- t(apply(fp.data.array,length(fp.dim),c))
  #eliminates the NA time-sequence from the dataset
  fp.idx.na           <- which(apply(fp.data.array,2,function(fp.x) any(is.na(fp.x))))
  fp.idxtovisit       <- 1:(fp.nx*fp.ny)
  if(length(fp.idx.na)>0) fp.idxtovisit <- fp.idxtovisit[-fp.idx.na]
  fp.data.array       <- fp.data.array[,fp.idxtovisit]
  fp.data.coord       <- fp.data.coord[fp.idxtovisit,]
  ## transforms the matrix into a matrix of the projections
  fp.data.array       <- RetProj()  

  ############################# LAUNCH ANALYSIS #############################  
  #for each pixel x time to visit do
  fp.lidxtovisit  <- length(fp.idxtovisit)
  if(fp.nproc==1){
    fp.res.visited  <- lapply(1:fp.lidxtovisit,function(fp.idx.pix){
      fp.neighb <- findNeighboors(fp.idx.pix)
      fp.neighb <- returnClosests(fp.idx.pix,fp.neighb)
      buildEstimate(fp.idx.pix,fp.neighb)
    })
  }else{
    fp.cl           <- makeCluster(getOption("cl.cores",fp.nproc))
    clusterExport(fp.cl, varlist=list("MultiTestH0"))  
    fp.res.visited  <- parLapply(fp.cl,1:fp.lidxtovisit,function(fp.idx.pix){
      fp.neighb <- findNeighboors(fp.idx.pix)
      fp.neighb <- returnClosests(fp.idx.pix,fp.neighb)
      buildEstimate(fp.idx.pix,fp.neighb)
    })
    stopCluster(fp.cl)}
    
  ############################# DATA REFORMATING #############################  
  #creates a list of size number of pixels x time, where the results will be stored
  fp.res.list      <- vector("list",fp.nx*fp.ny)
  for(fp.iterator in 1:fp.lidxtovisit){
    fp.idx.pix                                      <- fp.idxtovisit[fp.iterator]
    #creates the function fp.fctDenoise and initialises the parameters 
    fp.res.list[[fp.idx.pix]]                       <- fp.res.visited[[fp.iterator]]
    }
   
  return(fp.res.list)
  ### returns a list containing:
  ### - 'Vx' a vector containing all the indexes of the denoised pixel's neighboors
  ### - 'Ix' a vector containing the denoised signal
}, ex = function(){
# ## Not run:
# library(DynClust)
# data("adu340_4small",package="DynClust")
# #gain of the CCD camera
# #necessary in order to compute the variance of the dataset
# #estimated on calibration experiments
# G             <- 0.146
# #readout variance
# sro2          <- (16.4)^2
# #dataset's variance
# FT            <- adu340_4small
# FT_varhat     <- G*FT+G^2*sro2
# FT            <- FT/sqrt(FT_varhat)
# #launches the denoising step on the dataset with a statistical level of 5%
# denoisres     <- callDenoiseVoxel(FT,1,fp.mask.size=NULL,fp.nproc=2)
# #computes the average over time for each pixelxtime before denoising
# avg_before    <- apply(FT,1:2,mean)
# #computes the average over time for each pixelxtime after denoising
# avg_after     <- apply(getDenoisedSet(FT,denoisres),1:2,mean)
# #plotting results
# image(1:nrow(avg_before),1:ncol(avg_before),avg_before)
# x11();image(1:nrow(avg_after),1:ncol(avg_after),avg_after)
# ## End(Not run)
})
