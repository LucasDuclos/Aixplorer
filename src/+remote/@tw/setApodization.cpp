// #include <stdlib.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

// Please use only floats in this mex

#define M_PI_F 3.141592653589793238462643f
using namespace std;

/* ========================================================================== */
/* ========================================================================== */

// Function input arguments
#define	WFINI     prhs[1]
#define	APODFCT   prhs[2]
#define	DUTYCYCLE prhs[3]
#define	TXELEMTS  prhs[4]

// Function output arguments
#define	WFAPOD plhs[0]

/* ========================================================================== */
/* ========================================================================== */

int ControlSyntax(int /* nlhs */, int nrhs, const mxArray *prhs[])
{
    // Function variables
    std::string ErrStr("");
    size_t TxEldims, DutyCycledims;

    // There should be 4 input arguments
    if ( nrhs != 5 )
    {
        ErrStr += "The setApodization method";
        ErrStr += " of REMOTE.TW class needs 4 input";
        mexPrintf("%s\n", ErrStr.c_str());
        mexPrintf("The setApodization method of REMOTE.TW class needs 4 input arguments.\n");
        return 0;
    }

    // WFINI should belong to the single class
    if ( !mxIsClass(WFINI, "single") )
    {
        mexPrintf("Waveform should belong to the single class, and not %s.\n",
                mxGetClassName(WFINI));
        return 0;
    }

    // DUTYCYCLE should belong to the single class
    if ( !mxIsClass(DUTYCYCLE, "single") )
    {
        mexPrintf("Duty cycle should belong to the single class, and not %s.\n",
                mxGetClassName(DUTYCYCLE));
        return 0;
    }

    // TXELEMTS should belong to the int32 class
    if ( !mxIsClass(TXELEMTS, "int32") )
    {
        mexPrintf("TxElemts should belong to the single class, and not %s.\n",
                mxGetClassName(TXELEMTS));
        return 0;
    }

    // DutyCycle and TxElemts should have the same dimension
    DutyCycledims = (size_t) mxGetNumberOfElements(DUTYCYCLE);
    TxEldims      = (size_t) mxGetNumberOfElements(TXELEMTS);
    if ( DutyCycledims != TxEldims )
    {
        mexPrintf("TxElemts (%i) and DutyCycle (%i) should have the same number of elements.\n",
                DutyCycledims, TxEldims);
        return 0;
    }

    return 1;
}

/* ========================================================================== */
/* ========================================================================== */

float* buildApodWindow(int ApodFct, int FirstElem, int LastElem)
{
    // Function variables
    float *ApodWindow, *ApodPos;
    int NbElem ;
    
    // Initialize variables
    NbElem     = LastElem - FirstElem + 1;
    ApodPos    = new float[NbElem];
    ApodWindow = new float[NbElem];
    
    // Create the apodization window position
    for ( int k = 0; k < NbElem; k++ )
    {
        ApodPos[k] = 2 * (k - ((float) NbElem) / 2) / ((float) NbElem);
    }
    
    // Create the apodization window
    for ( int k = 0; k < NbElem; k++ )
    {
        if ( ApodFct == 1 )      // Bartlett
        {
            ApodWindow[k] = 1.0f - fabsf(ApodPos[k]);
        }
        else if ( ApodFct == 2 ) // Blackmanf
        {
            ApodWindow[k] = 21.0f/50.0f + 1.0f/2.0f * cosf(M_PI_F * ApodPos[k])
            + 2.0f/25.0f * cosf(2.0f * M_PI_F * ApodPos[k]);
        }
        else if ( ApodFct == 3 ) // Connes
        {
            ApodWindow[k] = powf(1.0f - powf(ApodPos[k], 2.0f), 2.0f);
        }
        else if ( ApodFct == 4 ) // Cosine
        {
            ApodWindow[k] = cosf(M_PI_F * ApodPos[k] / 2.0f);
        }
        else if ( ApodFct == 5 ) // Gaussian
        {
            ApodWindow[k] = expf( -powf(ApodPos[k], 2.0f));
        }
        else if ( ApodFct == 6 ) // Hamming
        {
            ApodWindow[k] = 27.0f/50.0f + 23.0f/50.0f * cosf(M_PI_F * ApodPos[k]);
        }
        else if ( ApodFct == 7 ) // Hanning
        {
            ApodWindow[k] = powf(cosf(M_PI_F * ApodPos[k] / 2.0f), 2.0f);
        }
        else if ( ApodFct == 8 ) // Welch
        {
            ApodWindow[k] = 1.0f - powf(ApodPos[k], 2.0f);
        }
        else
        {
            ApodWindow[k] = 1.0f;
        }
    }
    
    delete [] ApodPos;
    
    return ApodWindow;
}

/* ========================================================================== */
/* ========================================================================== */

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )
{
    
    // Function variables
    float *WFini, *WFapod, *ApodWindow, *DutyCycle;
    float TmpValue, Energy, Sign;
    mwSize *WFdims;
    size_t TxEldims;
    int *TxElemts;
    int ApodFct, kLoop, mLoop, nLoop, Nzeros1, Idx0;
    int FirstElem = 1000000000, LastElem = 0, TmpFirst = 0;
    
    /* ========================================================================== */
    
    // Control the syntax of the mex call
    if ( !ControlSyntax(nlhs, nrhs, prhs) )
    {
        mexErrMsgIdAndTxt("MATLAB:SETAPODIZATION",
                "The prerequisites of the method are not respected.");
    }

    // Matrixes dimensions
    WFdims   = (mwSize*) mxGetDimensions(WFINI);
    TxEldims = (size_t) mxGetNumberOfElements(TXELEMTS);

    // Create the matrixes
    WFAPOD = mxCreateNumericArray(2, WFdims, mxSINGLE_CLASS, (mxComplexity) !mxCOMPLEX);
    
    // Get data pointers
    WFini     = (float*) mxGetPr(WFINI);
    WFapod    = (float*) mxGetPr(WFAPOD);
    DutyCycle = (float*) mxGetPr(DUTYCYCLE);
    TxElemts  = (int*) mxGetPr(TXELEMTS);

    // Get apodization function
    ApodFct   = (int) mxGetScalar(APODFCT);

    // Estimate FirstElem and LastElem
    FirstElem = TxElemts[0] - 1;
    LastElem  = TxElemts[TxEldims-1] - 1;

    // Build the apodization aperture
    ApodWindow = buildApodWindow(ApodFct, FirstElem, LastElem);

    // Adapt the apodization window (duty cycle + firing elements)
    for ( kLoop = 0; kLoop < (int)TxEldims; kLoop++ )
    {
        Energy = 0.0f; Sign = 0.0f; TmpFirst = 0;
        for ( mLoop = 0; mLoop < WFdims[1]; mLoop++ )
        {
            TmpValue = WFini[(TxElemts[kLoop] - 1) + mLoop * WFdims[0]];
            if ( TmpValue != 0 )
            {
                if ( Sign == 0 )
                {
                    Sign     = TmpValue > 0 ? 1.0f : (TmpValue < 0.0f ? -1.0f : 0.0f);
                    Energy   = 0.0f;
                    TmpFirst = mLoop;
                }
                else if ( (TmpValue * Sign) < 0 )
                {
                    Energy  = floor(DutyCycle[kLoop] * Energy + 0.5f);
                    Nzeros1 = (int) ( (mLoop - TmpFirst - Energy) / 2.0f );
                    Idx0    = (TxElemts[kLoop] - 1) + (TmpFirst + Nzeros1) * WFdims[0];

                    for ( nLoop = 0; nLoop < Energy; nLoop++ ) // use <= if you need no 0 in the waveform, when not apodizing ?
                    {
                        WFapod[Idx0 + nLoop * WFdims[0]] = Sign;
                    }
                    Sign    = TmpValue > 0.0f ? 1.0f : (TmpValue < 0.0f ? -1.0f : 0.0f);
                    TmpFirst = mLoop;
                    Energy = 0.0f;
                }
                else
                {
                    Energy += ApodWindow[(TxElemts[kLoop] - 1 - FirstElem)]
                            * pow(WFini[(TxElemts[kLoop] - 1) + (mLoop-1) * WFdims[0]]
                            + WFini[(TxElemts[kLoop] - 1) + mLoop * WFdims[0]], 2) / 4;
                }
            }
            else if ( Sign != 0 )
            {
                Energy  = floor(DutyCycle[kLoop] * Energy + 0.5f);
                Nzeros1 = (int) ( (mLoop - TmpFirst - Energy) / 2.0f );
                Idx0    = (TxElemts[kLoop] - 1) + (TmpFirst + Nzeros1) * WFdims[0];

                for ( nLoop = 0; nLoop < Energy; nLoop++ ) {  // use <= if you need no 0 in the waveform, when not apodizing ?
                    WFapod[Idx0 + nLoop * WFdims[0]] = Sign;
                }
                Sign    = TmpValue > 0.0f ? 1.0f : (TmpValue < 0.0f ? -1.0f : 0.0f);
                TmpFirst = mLoop;
                Energy = 0.0f;
            }
        }
        
        if ( Energy > 0 )
        {
            Energy  = floor(DutyCycle[kLoop] * Energy + 0.5f);
            Nzeros1 = (int) ( (mLoop - TmpFirst - Energy) / 2.0f );
            Idx0    = (TxElemts[kLoop] - 1) + (TmpFirst + Nzeros1) * WFdims[0];

            for ( nLoop = 0; nLoop < Energy; nLoop++ ) {  // use <= if you need no 0 in the waveform, when not apodizing ?
                WFapod[Idx0 + nLoop * WFdims[0]] = Sign;
            }
        }
    }

    delete[] ApodWindow;
    
    return;
}
