#include "mex.h"
#include "complut.h"

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{

  /*  Check for proper number of arguments. */
  if (nrhs != 7)
      mexErrMsgTxt("Seven right parameters are required(pProbe, pAcqInfo, pReconInfo, pReconMux, pReconAbsc, pReconDelay, pOutputLut)");
  if (nlhs < 2)
      mexErrMsgTxt("2 or 3 left parameter is required (status).");


  /* Get pointers to the structure array fields and then get data. */

  double*   pProbe      = mxGetPr(prhs[0]) ;
  double*   pAcqInfo    = mxGetPr(prhs[1]) ;
  double*   pReconInfo  = mxGetPr(prhs[2]) ;
  int*      pReconMux   = (int *) mxGetData(prhs[3]);
  int*      pReconAbsc  = (int *) mxGetData(prhs[4]);
  double*   pReconDelay = mxGetPr(prhs[5]) ;
  

  // Set output pointer
  plhs[0]= mxCreateDoubleScalar(0);
  double* pStatus = mxGetPr(plhs[0]);

  int size_int = 0 ;
  /* call function to compute LUT size */
  CompLut(pProbe, pAcqInfo, pReconInfo, pReconMux, pReconAbsc, pReconDelay, (short *)&size_int);
  pReconInfo[19] = size_int ; // in bytes
  size_int /= 2 ;
  
  plhs[1] = mxCreateDoubleScalar(0);
  double* size = mxGetPr(plhs[1]);
  *size = size_int * 2 ;

  if (nlhs == 3)
  {
      mwSize ndim = 2 ;
      const mwSize dims[2] = { 1, size_int } ;
      plhs[2] = mxCreateNumericArray( ndim, dims, mxINT16_CLASS, mxREAL ) ;
  
      short*    pOutputLut  = (short *) mxGetData(plhs[2]);
  
      /* call function to compute LUT */
      pStatus[0] = CompLut(pProbe, pAcqInfo, pReconInfo, pReconMux, pReconAbsc, pReconDelay, pOutputLut);
  }
}

