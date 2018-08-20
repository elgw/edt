#include "mex.h"
#include "eudist2.c"


// Provides a MATLAB interface from eudist2.c

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[]) 
{

  if (nrhs<1 || nrhs>3) {
    mexErrMsgTxt("There should be 1 to three inputs.");
  }

  /* Check data types of the input arguments. */
  if (!(mxIsDouble(prhs[0]))) {
    mexErrMsgTxt("First argument must be of type double.");
  }

   double * V = (double *) mxGetPr(prhs[0]);
  mwSize VnDim = mxGetNumberOfDimensions(prhs[0]);
  const mwSize * Vdim = mxGetDimensions(prhs[0]);


  if(! (VnDim == 2))
    mexErrMsgTxt("Volume has to be 2D");

  plhs[0] = mxCreateNumericArray(2,  Vdim, mxDOUBLE_CLASS, mxREAL);

  double * D = (double *) mxGetPr(plhs[0]);

  edt2(V, D, Vdim[0], Vdim[1]);
 

}


