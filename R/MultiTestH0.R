MultiTestH0 <-
function
### 'MultiTestH0' tests if the vector or matrix 'fp.proj.matrix' given as an argument is significantly different from a null vector or matrix at a level alpha.
(
fp.proj.matrix,   
### a matrix of signals to be tested with as many columns as signals to test
fp.var,     
### a numeric indicating the variance of the dataset
fp.alpha=.05,       
### a numeric value indicating the level of the statistical multitest H0
fp.proc="bonferroni"         
### a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
### "fdr" method is not implemented yet
){
  #if fp.proj.matrix is a vector, then it will be reformated as a matrix with one column
  if(is.vector(fp.proj.matrix)) fp.proj.matrix  <- matrix(fp.proj.matrix,ncol=1)
  if(is.matrix(fp.var)) fp.proj.matrix  <- matrix(fp.proj.matrix,ncol=1)
  fp.proj.matrix[which(!is.finite(fp.proj.matrix))]  <- 0   
  #then normalize the data using this formula
  fp.proj.matrix  <- t(t(fp.proj.matrix)/sqrt(fp.var))
  fp.NCOL         <- ncol(fp.proj.matrix)
  fp.NROW         <- nrow(fp.proj.matrix)
  #sets the number of tests
  fp.NTEST        <- log2(fp.NROW+2)-1
  fp.proj.matrix  <-
.C("MultiTestH0_ccall",as.double(as.vector(fp.proj.matrix)),as.integer(fp.NROW),as.integer(fp.NCOL),as.integer(fp.NTEST),as.double(fp.alpha),as.integer(rep(0,ncol(fp.proj.matrix))),PACKAGE="DynClust")
  fp.proj.matrix  <- which(fp.proj.matrix[[6]]==0)
  return(fp.proj.matrix)
  ### returns TRUE if the test was not significant, FALSE if H1 is true (the mean is centered on 0) at a level alpha
  ### when 'fp.proj.matrix' is a matrix the function returns a vector of boolean values. 
}
