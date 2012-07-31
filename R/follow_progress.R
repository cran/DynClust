follow_progress <-
structure(function
### 'follow_progress' displays the progression of an iterative process
### and returns fp.prog a vector of values ranging from 0 to 100
### each value corresponding to a progression bar left for display 
(fp.idx.curr=0,
### a numeric indicating which is the current iteration 
fp.max.iter=0,
### a numeric indicating the total number of iterations
fp.prog=NULL,
### a vector of values ranging between 0 to 100 with a step of fp.step
### where each value corresponds to a progression bar to be displayed
### per default (when fp.prog=NULL) fp.prog = seq(0,100,by=fp.step)
fp.step=1 ##<< a numeric indicating the step necessary for constructing the fp.prog vector
##keyword<< internal  
){
  if(is.null(fp.prog)){
    fp.prog <- seq(0,100,by=fp.step) 
    cat(paste(c("0%",rep("|",length(fp.prog)-2),"100%","\n"),sep=""),sep="")
    return(fp.prog)}
  if(!all(range(fp.prog)>=0 & range(fp.prog)<=100)) stop("fp.prog is either NULL or a vector of values ranging between 0 and 100.")
  fp.val      <- trunc(trunc(100*fp.idx.curr/fp.max.iter)/fp.step)*fp.step
  fp.which.eq <- which(fp.val==fp.prog)
  if(length(fp.which.eq)==1){
    if(fp.prog[fp.which.eq]==0 |  fp.prog[fp.which.eq]==100){
      cat(paste(fp.prog[fp.which.eq],"%",sep=""))
      if(sum(fp.prog)!=sum(seq(0,100,by=fp.step))) cat("\n")
    }else{
      cat("|")
    }
    return(fp.prog[-fp.which.eq])
  }
  return(fp.prog)
  ### returns an updated version of fp.prog
}, ex = function(){
# ratiobars  <-  follow_progress()
#  for(i in 1:5000) ratiobars  <-  follow_progress(i,5000,ratiobars) 
})
