///////////////////////////////////////////////////////
//                                                   //
//  Copyright © 2005-2011 by SuperSonic Imagine, SA  //
//         Confidential - All Right Reserved         //
//                                                   //
///////////////////////////////////////////////////////



typedef struct float2 {
float x,y ;
} float2 ;

#include <stdio.h>
#include "beamformer.h"
#include "mex.h"
#include <time.h>

inline int min (int a,int b)
{
    if(a<b)
        return a ;
    else
        return b ;
}

int getIntParam(const mxArray * pIn, const char * fieldname, int & out)
{
    mxArray * pArray = NULL;
    pArray = mxGetField( pIn, 0, fieldname);
    if(pArray)
    {
        out = (int)mxGetScalar(pArray);
        mexPrintf("%s = %d\n",fieldname,out);
        return 0;
    }
    else
    {
        mexPrintf(fieldname);
        mexPrintf("not read \n");
        
        out = 0;
        return 1;
    }
}
int getFloatParam(const mxArray * pIn, const char * fieldname, float & out)
{
    mxArray * pArray = NULL;
    pArray = mxGetField( pIn, 0, fieldname);
    if(pArray)
    {
        out = (float)mxGetScalar(pArray);
        mexPrintf("%s = %g\n",fieldname,out);
        return 0;
    }
    else
    {
        mexPrintf(fieldname);
        mexPrintf("not read \n");
        out = 0;
        return 1;
    }
}

bool copyBeamformingStruct(const mxArray * pIn, SBeamformParameters * out)
{

    // float parameters

    if(getFloatParam(pIn,"speedOfSound",out->m_speedOfSound))
        return false;
    
    out->m_invSpeedOfSound = 1.0f/out->m_speedOfSound;
    
    
    if(getFloatParam(pIn,"ProbeRadius",out->m_ProbeRadius))
        return false;
    if(getFloatParam(pIn,"piezoPitch",out->m_piezoPitch))
        return false;
    if(getFloatParam(pIn,"demodFreq",out->m_demodFreq))
        return false;
    if(getFloatParam(pIn,"rOrigin",out->m_rOrigin))
        return false;
    if(getFloatParam(pIn,"fNumber",out->m_fNumber))
        return false;
    if(getFloatParam(pIn,"peakDelay",out->m_peakDelay))
        return false;
    if(getFloatParam(pIn,"reconPitch",out->m_reconPitch))
        return false;
    if(getFloatParam(pIn,"linePitch",out->m_linePitch))
        return false;
    if(getFloatParam(pIn,"pixelPitch",out->m_pixelPitch))
        return false;
    if(getFloatParam(pIn,"lambda",out->m_lambda))
        return false;
    if(getFloatParam(pIn,"xOrigin",out->m_xOrigin))
        return false;

    // integer parameters
    if(getIntParam(pIn,"nbPiezos",out->m_nbPiezos))
        return false;
    if(getIntParam(pIn,"nbAngles",out->m_nbAngles))
        return false;
    if(getIntParam(pIn,"channelOffset",out->m_channelOffset))
        return false;
    if(getIntParam(pIn,"firstSample",out->m_firstSample))
        return false;
    if(getIntParam(pIn,"nbSamples",out->m_nbSamples))
        return false;
    if(getIntParam(pIn,"nbLinesPerRecon",out->m_nbLinesPerRecon))
        return false;
    if(getIntParam(pIn,"nbPixelsPerLine",out->m_nbPixelsPerLine))
        return false;
    if(getIntParam(pIn,"nbRecon",out->m_nbRecon))
        return false;
    if(getIntParam(pIn,"normMode",out->m_normMode))
        return false;
    if(getIntParam(pIn,"synthAcq",out->m_synthAcq))
        return false;
    if(getIntParam(pIn,"firstAcquiredChannel",out->m_firstAcquiredChannel))
        return false ;
    if(getIntParam(pIn,"firstRecon",out->m_firstRecon))
        return false ;
    if(getIntParam(pIn,"usegpu",out->m_usegpu))
        return false ;
    
   // angles

    mxArray * pArray;
    pArray = mxGetField( pIn, 0, "angles");
    if(pArray)       
    {
        double * pDouble = mxGetPr(pArray);
        for(int i=0; i< out->m_nbAngles;i++)
        {
            out->m_angles[i] = (float)pDouble[i];
            mexPrintf("out->m_angles[i] =%g\n",out->m_angles[i]);
        }

    }
    else
    {
         return false;
    }
    
    // time origin
    pArray = mxGetField( pIn, 0, "timeOrigin");
    if(pArray)       
    {
        double * pDouble = mxGetPr(pArray);
        for(int i=0 ; i< out->m_nbAngles;i++)
        {
            out->m_timeOrigin[i] = (float)pDouble[i];
            mexPrintf("out->m_timeOrigin[i] =%g\n",out->m_timeOrigin[i]);
        }
        for(int i=out->m_nbAngles ; i< NB_MAX_ANGLES;i++)
        {
            out->m_timeOrigin[i] = 0.0f;
        }
    }
    else
    {
         return false ;
    }
    
    return true;
}

void printBeamformingStruct(SBeamformParameters * out)
{

    // float parameters

    if(getFloatParam(pIn,"speedOfSound",out->m_speedOfSound))
    
    out->m_invSpeedOfSound = 1.0f/out->m_speedOfSound;
    
    
    mexPrintf("ProbeRadius=%g", out->m_ProbeRadius))
    mexPrintf("piezoPitch=%g", out->m_piezoPitch))
    mexPrintf("demodFreq=%g", out->m_demodFreq))
    mexPrintf("rOrigin=%g", out->m_rOrigin))
    mexPrintf("fNumber=%g", out->m_fNumber))
    mexPrintf("peakDelay=%g", out->m_peakDelay))
    mexPrintf("reconPitch=%g", out->m_reconPitch))
    mexPrintf("linePitch=%g", out->m_linePitch))
    mexPrintf("pixelPitch=%g", out->m_pixelPitch))
    mexPrintf("lambda=%g", out->m_lambda))
    mexPrintf("xOrigin=%g", out->m_xOrigin))

    // integer parameters
    mexPrintf("nbPiezos=%d", out->m_nbPiezos))
    mexPrintf("nbAngles=%d", out->m_nbAngles))
    mexPrintf("channelOffset=%d", out->m_channelOffset))
    mexPrintf("firstSample=%d", out->m_firstSample))
    mexPrintf("nbSamples=%d", out->m_nbSamples))
    mexPrintf("nbLinesPerRecon=%d", out->m_nbLinesPerRecon))
    mexPrintf("nbPixelsPerLine=%d", out->m_nbPixelsPerLine))
    mexPrintf("nbRecon=%d", out->m_nbRecon))
    mexPrintf("normMode=%d", out->m_normMode))
    mexPrintf("synthAcq=%d", out->m_synthAcq))
    mexPrintf("firstAcquiredChannel=%d",  out->m_firstAcquiredChannel))
    mexPrintf("firstRecon=%d", out->m_firstRecon))
    mexPrintf("usegpu=%d", out->m_usegpu))
    
    // angles
    for(int i=0; i< out->m_nbAngles;i++)
    {
        mexPrintf("out->m_angles[i] =%g\n",out->m_angles[i]);
    }
    
    // time origin
    for(int i=0 ; i< out->m_nbAngles;i++)
    {
        mexPrintf("out->m_timeOrigin[i] =%g\n",out->m_timeOrigin[i]);
    }
    
}

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{

    time_t start,end;
    bool status = 0;
    SBeamformParameters bfStruct;

    // check the arguments
    if(nrhs != 3)
        mexErrMsgTxt("3 input required");

    const mxArray * pInStruct = prhs[0];
    const mxArray * pRFArray = prhs[1];
    const mxArray * palignedOffset = prhs[2];
    
    const int alignedOffset = ((int)mxGetScalar(palignedOffset))/2;
    const int bufferSize = mxGetNumberOfElements(pRFArray)-alignedOffset;
    short * pInRF = (short*)mxGetData(pRFArray)+alignedOffset;
    
    // parse structure
    //mexPrintf("alignedOffset = %d\n",alignedOffset);
    status = copyBeamformingStruct(pInStruct, &bfStruct);
    printBeamformingStruct(&bfStruct);
    if(status)
    {
        const int nbEffectiveChannels = min(bfStruct.m_nbPiezos,128);
        const int reducedBufferSize = bfStruct.m_channelOffset*nbEffectiveChannels;
        // create output matrix
        plhs[0] = mxCreateNumericMatrix(2*bfStruct.m_nbPixelsPerLine, bfStruct.m_nbRecon*bfStruct.m_nbLinesPerRecon, mxSINGLE_CLASS, mxREAL);
        float2 * pOutImageIQ = (float2*)mxGetData(plhs[0]);      
        float elapsedTimeMS ;
        if(bfStruct.m_usegpu)
        {
            start = clock();
            cudaIQBeamform(&bfStruct,  pInRF , reducedBufferSize , pOutImageIQ,elapsedTimeMS);
            end = clock();
            double d =  ((double) (end - start)) / CLOCKS_PER_SEC;
            mexPrintf ("global call of beamformer: cpu timeu is  %f ms \n", d*1e3);
        }
        else
            iq_beamform_kernel_cpu(pInRF, bufferSize, pOutImageIQ, bfStruct);
        
        mexPrintf("kernel call (without cuda init an gpu malloc) is %f ms\n",elapsedTimeMS);
    }
    else {
        mexErrMsgTxt("problem with a parameter");
    }
}
