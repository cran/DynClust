getChildren <-
function
### 'getChildren' returns the indexes of all the direct and indirect 
### neighboors of the pixel x time at index fp.idx.cur.pix
### which are still inside the list of pixels x time left to clust.
(fp.idx.cur.pix,
### a numeric indicating the index of the current pixel x time
fp.pixtoclust,
### a vector listing the pixels x time left to clust 
fp.res.denois
### a list containing the results of the denoising procedure (cf. callDenoiseVoxel)
##keyword<< internal
){
  #sets N as Vx, the neighboorhood of x
  fp.N    <- fp.res.denois[[fp.idx.cur.pix]]$Vx
  #sets an index to 0
  fp.idx  <- 0
  #while this index is smaller than the size of fp.N do
  while(fp.idx<length(fp.N)){
    #increments the index
    fp.idx          <- fp.idx+1
    #merges the N to the neighboorhood Vx of the pixel x time at index fp.idx (+elimination of repeated values)
    fp.Ntemp        <- unique(c(fp.N,fp.res.denois[[fp.idx]]$Vx))
    #checks for the neighboors that are still in the list of pixels x time to clust
    fp.neighbmatch  <- which(fp.Ntemp %in% fp.pixtoclust)
    #selects only the available neighboors 
    fp.N            <- fp.Ntemp[fp.neighbmatch]
}
#returns fp.N
fp.N
### returns a vector containing the indexes of all the direct and indirect neighboors of the pixel x time at index fp.idx.cur.pix
### which are still inside the list of pixels x time left to clust
}
