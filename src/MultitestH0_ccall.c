#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

/**
 * MultiTestH0_ccall(double *MAT, int *NROW, int *NCOL,int *NBTEST,double *ALPHA,int *RES)
 *
 * @brief This function applies a multitest H0 
 *
 * @param *MAT a pointer to the  in R.
 * @param *NROW a pointer to the  in R.
 * @param *NCOL a pointer to the  in R.
 * @param *NTEST a pointer to .
 * @param *ALPHA a pointer to .
 * @param *winners a pointer to .
 *
 */


int Winner(double *tab,int size,double thres){
  double sum = 0;
  for(int irow=0; irow<=(size-1);irow=irow+1){
    sum = sum+tab[irow]*tab[irow];
    if(2*sum>thres) return(1);
  }    
  return(0);
}
void MultiTestH0_ccall(double *MAT, int *NROW, int *NCOL, int *NTEST,double *ALPHA,int *winners){
  int ntest = *NTEST;
  int jcol = 0;
  double *tab;
  double thresh=0;
  int size=2^ntest;
  int sum = 0;
  int starting = 0;
  while(ntest>=1){
    size    = R_pow(2,ntest);
    thresh  = qchisq(1-((*ALPHA)/((*NTEST)+1)),size,1,0);
    for(jcol=0;jcol<=(*NCOL-1);jcol++){
      if(winners[jcol]==0){
        tab             = &MAT[jcol*(*NROW)+starting];
        winners[jcol]   = Winner(tab,size,thresh); 
      }
      sum = sum+winners[jcol];
      if(sum==*NCOL) break;
    }
    sum      = 0;
    starting = starting+size;
    ntest    = ntest-1;
  }

}

