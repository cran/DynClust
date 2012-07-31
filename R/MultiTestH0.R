MultiTestH0 <-
function
### 'MultiTestH0' tests if the vector or matrix 'fp.diff.matrix' given as an argument is significantly different from a null vector or matrix at a level alpha.
(
fp.diff.matrix,   
### a matrix of signals to be tested with as many columns as signals to test
fp.stab.var,     
### a numeric or array indicating the variance of the dataset
fp.alpha=.05,       
### a numeric value indicating the level of the statistical multitest H0
fp.proc="bonferroni"         
### a character either "bonferroni" or "fdr" indicating which method to use for the multitest H0
### "fdr" method is not implemented yet
){
  #if fp.diff.matrix is a vector, then it will be reformated as a matrix with one column
  if(is.vector(fp.diff.matrix)) fp.diff.matrix  <- matrix(fp.diff.matrix,ncol=1)
  fp.stab.var[which(!is.finite(fp.diff.matrix))]     <- 0
  fp.diff.matrix[which(!is.finite(fp.diff.matrix))]  <- 0
  fp.stab.var[which(!is.finite(fp.stab.var))]        <- 0   
  #sets the number of tests
  fp.ntest    <- log2(nrow(fp.diff.matrix))
  #compare the number of rows (length of each signal) with the number of tests 
  if(fp.ntest!=round(fp.ntest)){stop("the time dimension is not of the form n=2^d => log2(n) is not an integer")}
  #then normalize the data using this formula
  fp.diff.matrix  <- fp.diff.matrix/sqrt(fp.stab.var) 
  res <- .C("MultiTestH0_ccall",as.double(as.vector(fp.diff.matrix)),as.integer(nrow(fp.diff.matrix)),as.integer(ncol(fp.diff.matrix)),as.integer(fp.ntest),as.double(fp.alpha),as.integer(rep(1,ncol(fp.diff.matrix))),PACKAGE="DynClust")
  ret1  <- which(res[[6]]==1)
  return(ret1)
  ### returns TRUE if the test was not significant, FALSE if H1 is true (the mean is centered on 0) at a level alpha
  ### when 'fp.diff.matrix' is a matrix the function returns a vector of boolean values. 
}
