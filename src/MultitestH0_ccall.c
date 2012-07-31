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
  int i = 0, j = 0,compteur=0, nbstep,taille = ((*NROW)*(*NCOL)), newnrow=(*NROW);
  double norm2 = 0, tresh = 0;
  while(newnrow >= 2){
    tresh   = qchisq(1-((*ALPHA)/(*NBTEST)),newnrow/2,1,0);
    j       = 0;
    for(i=0; i<taille-1; i=i+2){
      MAT[i/2]  = (MAT[i]+MAT[i+1])/2;
      norm2   = norm2+MAT[i/2]*MAT[i/2];
      if(i==(newnrow-2)+newnrow*j){
        norm2 = norm2*2;
        if(norm2>tresh) RES[j] = 0;
        j++;
        norm2=0; 
      }
    }
    taille  = taille/2;
    newnrow = newnrow/2;
  }
}

