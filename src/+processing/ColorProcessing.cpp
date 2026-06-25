#include <math.h>
#include "mex.h"
#include "matrix.h"

/* ========================================================================== */
/* ========================================================================== */

// Function input arguments
#define	IQ              prhs[0]
#define	MASK            prhs[1]
#define	R0REFVALUE      prhs[2]
#define	R0COMPRESSIONDB prhs[3]

// Function output arguments
#define	WFIQ   plhs[0]
#define	R0     plhs[1]
#define	R1     plhs[2]
#define	STATUS plhs[3]

/* ========================================================================== */
/* ========================================================================== */

int ControlSyntax(int nlhs, int nrhs, const mxArray *prhs[])
{
    
    // Function variables
    int IQnd, WFIQnd, MaskNd;
    mwSize *IQdims, *WfIQdims, *MaskDims;


/* ========================================================================== */

    // Control the number of arguments
    if ( (nlhs != 3) && (nlhs != 4) )
    {
        mexPrintf("There should be three output arguments.\n");
        return 0;
    }
    else if ( nrhs != 4 )
    {
        mexPrintf("Four input arguments are required.\n");
        return 0;
    }
    
/* ========================================================================== */
    
    // Control the IQ image
    if ( mxIsComplex(IQ) ) // IQ needs to be complex
    {
        IQnd = mxGetNumberOfDimensions(IQ);
        if ( IQnd != 3 ) // IQ needs to be 3D
        {
            mexPrintf("IQ data need to be a 3-D complex array, and not %d-D.\n", 
                IQnd);
            return 0;
        }
    }
    else
    {
        mexPrintf("IQ data need to be a 3-D complex array.\n");
        return 0;
    }
    IQdims = (mwSize*) mxGetDimensions(IQ);
    
    // IQ data should belong to the single class
    if ( !mxIsClass(IQ, "single") )
    {
        mexPrintf("IQ data need to belong to the single class, and not %s.\n",
                mxGetClassName(MASK));
        return 0;
    }
    
/* ========================================================================== */
    
    // Control the MASK
    MaskNd = mxGetNumberOfDimensions(MASK);
    if ( MaskNd != 2 ) // MASK needs to be 2D
    {
        mexPrintf("Mask data need to be a 2-D array, and not %d-D.\n", MaskNd);
        return 0;
    }
    
    // Control the dimensions of the mask
    MaskDims = (mwSize*) mxGetDimensions(MASK);
    if ( MaskDims[0] != MaskDims[1] ) // MASK needs to be a square matrix
    {
        mexPrintf("The Mask data should be a squared matrix and not %dx%d.\n",
            MaskDims[0], MaskDims[1]);
        return 0;
    }
    else
    {
        bool Test = ceil((float) IQdims[2] / MaskDims[0]) 
                        == floor((float) IQdims[2] / MaskDims[0]);
        if ( !Test )
        {
            mexPrintf(
                "The IQ 3rd dimension (%d) is not a Mask size multiple (%d).\n",
                IQdims[2], MaskDims[0]);
            return 0;
        }
    }
    
    // MASK data should belong to the single class
    if ( !mxIsClass(MASK, "single") )
    {
        mexPrintf("MASK data need to belong to the single class, and not %s.\n",
                mxGetClassName(MASK));
        return 0;
    }
    
/* ========================================================================== */

    
    
/* ========================================================================== */
    
    return 1;
    
}

/* ========================================================================== */
/* ========================================================================== */

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )
{
    
    // Function variables
    float *IQr, *IQi, *Mask, *WfIQr, *WfIQi, *R0r, *R0i, *R1r, *R1i,
            R0RefValue, R0CompressiondB;
    int IQnElmts, RnElmts, k, m, IdM, IdM0, n, IdN;
    mwSize *IQdims, *MaskDims, Rdims[3];
    
/* ========================================================================== */
    
    // Control the syntax of the mex call
    if ( !ControlSyntax(nlhs, nrhs, prhs) )
    {
        mexErrMsgIdAndTxt("MATLAB:COLORPROCESSING", 
                "The prerequisites of the method are not respected.");
    }
    else // create the WFIQ, R0 and R1 matrixes
    {
        // Matrixes dimensions
        IQdims   = (mwSize*) mxGetDimensions(IQ);
        MaskDims = (mwSize*) mxGetDimensions(MASK);
        Rdims[0] = IQdims[0];
        Rdims[1] = IQdims[1];
        Rdims[2] = IQdims[2] / MaskDims[0];
        
        // Create the matrixes
        WFIQ = mxCreateNumericArray(3, IQdims, mxSINGLE_CLASS, mxCOMPLEX);
        R0   = mxCreateNumericArray(3, Rdims, mxSINGLE_CLASS, (mxComplexity) !mxCOMPLEX);
        R1   = mxCreateNumericArray(3, Rdims, mxSINGLE_CLASS, mxCOMPLEX);
    }
    
/* ========================================================================== */
    
    // Get data pointers
    Mask  = (float*) mxGetPr(MASK);
    IQr   = (float*) mxGetPr(IQ);   IQi   = (float*) mxGetPi(IQ);
    WfIQr = (float*) mxGetPr(WFIQ); WfIQi = (float*) mxGetPi(WFIQ);
    R0r   = (float*) mxGetPr(R0);
    R1r   = (float*) mxGetPr(R1);   R1i   = (float*) mxGetPi(R1);

    // Get display parameters
    R0RefValue      = mxGetScalar(R0REFVALUE);
    R0CompressiondB = mxGetScalar(R0COMPRESSIONDB);
    
/* ========================================================================== */
    
    // Loop over all pixels
    IQnElmts = IQdims[0] * IQdims[1] * IQdims[2];
    RnElmts  = Rdims[0] * Rdims[1] * Rdims[2];
    
    for ( k = 0; k < RnElmts; k++ )
    {
        // Loop over the WFIQ images
        for ( m = 0; m < MaskDims[0]; m++ )
        {
            // Estimate the WFIQ index
            IdM = fmod((float) k, (IQdims[0] * IQdims[1])) + IQdims[0] * IQdims[1]
                    * (m +  floor((float) k / (IQdims[0] * IQdims[1])) * MaskDims[0]);
            
            // Loop over the IQ images
            for ( n = 0; n < MaskDims[0]; n++ )
            {
                // Estimate the IQ index
                IdN = IdM + (n - m) * (IQdims[0] * IQdims[1]);
                
                // Update the WFIQ images
                WfIQr[IdM] += IQr[IdN] * Mask[m + n * MaskDims[0]];
                WfIQi[IdM] += IQi[IdN] * Mask[m + n * MaskDims[0]];
            }
            
            // Update the R0 images
            R0r[k] += WfIQr[IdM] * WfIQr[IdM] + WfIQi[IdM] * WfIQi[IdM];
            
            // Update the R1 images
            if ( m > 0 )
            {
                IdM0 = IdM - (IQdims[0] * IQdims[1]);
                R1r[k] += WfIQr[IdM] * WfIQr[IdM0] + WfIQi[IdM] * WfIQi[IdM0];
                R1i[k] += WfIQr[IdM] * WfIQi[IdM0] - WfIQi[IdM] * WfIQr[IdM0];
            }
        }
        
        // Normalize the R0 and R1 images
        if ( R0r[k] > (MaskDims[0] * R0RefValue) )
        {
            // Update R0
            if ( R0r[k] > ( (MaskDims[0] * R0RefValue)
                                            * pow(10, R0CompressiondB / 10) ) )
            {
                R0r[k] = R0CompressiondB;
            }
            else
            {
                R0r[k]  = 10 * log10(R0r[k] / (MaskDims[0] * R0RefValue));
            }
            
            // Update R1
            R1r[k] /= (MaskDims[0] - 1);
            R1i[k] /= (MaskDims[0] - 1);
        }
        else
        {
            // Update R0
            R0r[k] = 0;
            
            // Update R1
            R1r[k] = 0;
            R1i[k] = 0;
        }
    }
    
/* ========================================================================== */
/* ========================================================================== */

    return;
}
