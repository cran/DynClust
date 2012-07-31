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
 * @param *NBTEST a pointer to .
 * @param *ALPHA a pointer to .
 * @param *RES a pointer to .
 *
 */
void MultiTestH0_ccall(double *MAT, int *NROW, int *NCOL,int *NBTEST,double *ALPHA,int *RES){
  int i = 0, j = 0, taille = ((*NROW)*(*NCOL)), newnrow=(*NROW);
  double norm2 = 0, tresh = 0;
 //while the time dimension has at least 2 observations do  
  while(newnrow >= 2){
    tresh   = qchisq(1-((*ALPHA)/(*NBTEST)),newnrow/2,1,0);
    j       = 0;
    //for each pair of observation
    for(i=0; i<taille-1; i=i+2){
      if(RES[j] == 1){
        //compute the mean of the pair of observations
        MAT[i/2]  = (MAT[i]+MAT[i+1])/2;
        //compute the sum of the square means
        norm2   = norm2+MAT[i/2]*MAT[i/2];
        //if 2 times the sum of the square means is greater than the threshold increment i
        if(norm2*2>tresh) i=(newnrow-2)+newnrow*j;     
        //if i has reach the end of a virtual column
        if(i==(newnrow-2)+newnrow*j){
          if(norm2*2>tresh){
            RES[j] = 0;
          }
          j++;
          norm2=0;}
        
      }else{
        i=(newnrow-2)+newnrow*j; 
        j++;
        norm2=0;
      }
      
    }
    taille  = taille/2;
    newnrow = newnrow/2;
  }
}



