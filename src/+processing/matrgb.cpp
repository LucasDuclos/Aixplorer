#include <math.h>
#include "mex.h"
#include "matrix.h"

/* ========================================================================== */
/* ========================================================================== */

// Function input arguments
#define	IMMATRIX   prhs[0]
#define	CMAP       prhs[1]
#define	CLIM       prhs[2]

// Function output arguments
#define	IMRED      plhs[0]
#define	IMGREEN    plhs[1]
#define	IMBLUE     plhs[2]

/* ========================================================================== */
/* ========================================================================== */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    // Function variables
    float *ImMatrix, *CLim, tmp;
    unsigned char *ImRed, *ImGreen, *ImBlue, *CMap;
    mwSize mCmap, nCmap, *dims;
    int k, nElmts;
    
    // Control the number of arguments
    if (nrhs != 3)
    {
        mexErrMsgTxt("Three input arguments are required...");
    }
    else if (nlhs != 3)
    {
        mexErrMsgTxt("Three output arguments are required...");
    }
    
/* ========================================================================== */
    
    // Determine the output matrix dimensions
    nElmts = mxGetNumberOfElements(IMMATRIX);
    mCmap  = mxGetM(CMAP);
    nCmap  = mxGetN(CMAP);
    dims   = (mwSize*) mxGetDimensions(IMMATRIX);
    
    // Create matrix for output
    IMRED   = mxCreateNumericArray(2, dims, mxUINT8_CLASS, (mxComplexity) !mxCOMPLEX);
    IMGREEN = mxCreateNumericArray(2, dims, mxUINT8_CLASS, (mxComplexity) !mxCOMPLEX);
    IMBLUE  = mxCreateNumericArray(2, dims, mxUINT8_CLASS, (mxComplexity) !mxCOMPLEX);
    
    // Assign input pointers
    ImMatrix = (float*) mxGetPr(IMMATRIX);
    CMap     = (unsigned char*) mxGetPr(CMAP);
    CLim     = (float*) mxGetPr(CLIM);
    
    // Assign output pointers
    ImRed   = (unsigned char*) mxGetPr(IMRED);
    ImGreen = (unsigned char*) mxGetPr(IMGREEN);
    ImBlue  = (unsigned char*) mxGetPr(IMBLUE);
    
/* ========================================================================== */
    
    // Loop to evaluate the values for each channel
    for (k = 0; k < nElmts; k++)
    {
        
        // Format values
        tmp = ((ImMatrix[k] - CLim[0]) * ((float) mCmap-1) / (CLim[1] - CLim[0]));
        if (tmp > (mCmap-1)) { tmp = (mCmap-1); }
        if (tmp < 0) { tmp = 0; }
        
        // Format color values
        ImRed[k]   = CMap[(unsigned short int) tmp];
        ImGreen[k] = CMap[(unsigned short int) tmp + mCmap];
        ImBlue[k]  = CMap[(unsigned short int) tmp + 2 * mCmap];

    }
    
    return;
}
