///////////////////////////////////////////////////////
//                                                   //
//  Copyright © 2005-2011 by SuperSonic Imagine, SA  //
//         Confidential - All Right Reserved         //
//                                                   //
///////////////////////////////////////////////////////
//
// cuda/beamformer_kernel.cu
//

typedef struct float2 {
float x,y ;
} float2 ;

#include <stdio.h>
#include <math.h>
#include "beamformer.h"
#include "mex.h"

inline int min(int a, int b)
{
    if(a<b)
        return a;
    else
        return b;
}

inline int max(int a, int b)
{
    if(a>b)
        return a;
    else
        return b;
}


float apod_c(int n, int apertureLength)
{
    return (0.54f - 0.46f * cos(2.0f * (float)(M_PI) * (float)(n+1) / ((float)apertureLength - 1.0f)));
}

void iq_beamform_kernel_cpu(const short* pInRF,const int bufferSize, float2* pOutImageIQ, SBeamformParameters g_bf_params)
{
    // output pixel coords
    
    for(int ix=0; ix< g_bf_params.m_nbRecon; ix++)
    {
        float thetaFiring = g_bf_params.m_firstFiringAngle + ix*g_bf_params.m_firingAnglePitch;
        for(int ir=0; ir<g_bf_params.m_nbPixelsPerLine; ir++)
        {

		float RfromApex = (float)ir * g_bf_params.m_pixelPitch + g_bf_params.m_rOrigin;
        float aperture = RfromApex/g_bf_params.m_fNumber;
		int halfAperture = (int) round(0.5f * aperture / g_bf_params.m_piezoPitch);
        int firstChannel0 = g_bf_params.m_nbPiezos/2 - halfAperture; 
        int lastChannel0 = g_bf_params.m_nbPiezos/2 + halfAperture;
		
		//int columnOffset = ix * g_bf_params.m_nbLinesPerRecon * g_bf_params.m_nbPixelsPerLine;

		//float refDelay = 2.0f * Z * g_bf_params.m_invSpeedOfSound + g_bf_params.m_peakDelay;
		
		const float w0 = -1.0f/6.0f;
		const float w1 = 0.5f;
		const float w2 = -0.5f;
		const float w3 = 1.0f/6.0f;
		const float eps = 0.000001f;
		const float normal_fudge = 1.0f;

		int firstChannel = 0;
		int lastChannel  = min(NB_ACQ_CHANNELS,g_bf_params.m_nbPiezos);
				
		float normCoef = eps;
		float IQ[2];
		int offset;

		for (int channel = firstChannel; channel < lastChannel; channel++)
		{
	#ifdef USE_CHANNEL_OFFSET
		    offset = g_bf_chanOffsets[channel];
	#else
		    offset = (channel % NB_ACQ_CHANNELS) * g_bf_params.m_channelOffset + ix * g_bf_params.m_nbSamples;
		   
	#endif
		    float apodCoef = apod_c(channel - firstChannel0, lastChannel0 - firstChannel0+1);
		    float xPiezo = (float)(channel +0.5 - g_bf_params.m_nbPiezos/2) * g_bf_params.m_piezoPitch;
            float zPiezo = g_bf_params.m_zApex;
		    		    
            const short* channelStartPtr = pInRF + offset;
            
            #pragma unroll 2
            for (int iPixel = 0; iPixel < g_bf_params.m_nbLinesPerRecon; iPixel++)
            {
                int columnOffset = (ix*g_bf_params.m_nbLinesPerRecon+iPixel) * g_bf_params.m_nbPixelsPerLine;
                float thetaLine = thetaFiring + iPixel * g_bf_params.m_linePitch + g_bf_params.m_thetaOrigin;
                float RfromSkinLine = RfromApex - zPiezo/cos(thetaLine);
                float X = RfromApex*sin(thetaLine);
                float Z = RfromApex*cos(thetaLine);
                
                float forwardDelay = RfromSkinLine*g_bf_params.m_invSpeedOfSound;
                float returnDelay  = sqrt((Z-zPiezo) * (Z-zPiezo) + (X - xPiezo) * (X - xPiezo)) * g_bf_params.m_invSpeedOfSound;
                float timeDelay = forwardDelay + returnDelay + g_bf_params.m_timeOrigin[ix];

                int delay = (int) floor(timeDelay * g_bf_params.m_demodFreq) - g_bf_params.m_firstSample/2;
                
                if(delay < 0)
                {
                    mexPrintf("forwardDelay = %g, returnDelay = %g , ix = %d \n",forwardDelay,returnDelay,ix);
                    mexErrMsgTxt("error: delay < 0");
                }                    
                
                if( 2*delay-2+1 + offset >= bufferSize )
                {
                    mexPrintf("forwardDelay = %g, returnDelay = %g , ix = %d \n",forwardDelay,returnDelay,ix);
                    mexErrMsgTxt("error: delay > acquisition time");
                }
                float deltaDelay = 2.0f * (float)(M_PI) * (timeDelay) * g_bf_params.m_demodFreq;
                float phi = timeDelay * g_bf_params.m_demodFreq - floor(timeDelay * g_bf_params.m_demodFreq);
#ifdef USE_BRANCH
                if (phi < eps)
                {
                    IQ[0] = channelStartPtr[2*delay-2];
                    IQ[1] = channelStartPtr[2*delay-2+1];
                }
                else
                {
                    float X0 = phi+1.0f;
                    float X1 = phi;
                    float X2 = phi-1.0f;
                    float X3 = phi-2.0f;
                    float l = X0*X1*X2*X3;
                    IQ[0] = l*(w0*channelStartPtr[2*delay-2]/X0 + w1*channelStartPtr[2*delay]/X1 + w2*channelStartPtr[2*delay+2]/X2 + w3*channelStartPtr[2*delay+4]/X3);
                    IQ[1] = l*(w0*channelStartPtr[2*delay-2+1]/X0 + w1*channelStartPtr[2*delay+1]/X1 + w2*channelStartPtr[2*delay+2+1]/X2 + w3*channelStartPtr[2*delay+4+1]/X3);
                }
#else
                // If the previous if/else statement is too slow
                // You can replace it by the following code
                float X0 = phi+1.0f;
                float X1 = phi+eps;
                float X2 = phi-1.0f;
                float X3 = phi-2.0f;

                // expand some terms and eliminate divides
                float l0 = X1*X2*X3*w0;
                float l1 = X0*X2*X3*w1;
                float l2 = X0*X1*X3*w2;
                float l3 = X0*X1*X2*w3;
                IQ[0] = l0*channelStartPtr[2*delay-2] + l1*channelStartPtr[2*delay] + l2*channelStartPtr[2*delay+2] + l3*channelStartPtr[2*delay+4];
                IQ[1] = l0*channelStartPtr[2*delay-2+1] + l1*channelStartPtr[2*delay+1] + l2*channelStartPtr[2*delay+2+1] + l3*channelStartPtr[2*delay+4+1];
    #endif
//                 pOutImageIQ[columnOffset + ir + iPixel].x += (IQ[0]*cos(deltaDelay) + IQ[1]*sin(deltaDelay))*apodCoef;
//                 pOutImageIQ[columnOffset + ir + iPixel].y += (-IQ[0]*sin(deltaDelay) + IQ[1]*cos(deltaDelay))*apodCoef;
                pOutImageIQ[columnOffset + ir].x += (IQ[0]*cos(deltaDelay) + IQ[1]*sin(deltaDelay))*apodCoef;
                pOutImageIQ[columnOffset + ir].y += (-IQ[0]*sin(deltaDelay) + IQ[1]*cos(deltaDelay))*apodCoef;
            }

		    normCoef += apodCoef * apodCoef;
		}

		if (g_bf_params.m_normMode < 100)
		{
		    #pragma unroll 2
		    for (int iPixel = 0; iPixel < g_bf_params.m_nbLinesPerRecon; iPixel++)
		    {
                int columnOffset = (ix*g_bf_params.m_nbLinesPerRecon+iPixel) * g_bf_params.m_nbPixelsPerLine;
// 		        pOutImageIQ[columnOffset + ir + iPixel].x /= normCoef;
// 		        pOutImageIQ[columnOffset + ir + iPixel].y /= normCoef;
// 		        pOutImageIQ[columnOffset + ir + iPixel].x *= normal_fudge;
// 		        pOutImageIQ[columnOffset + ir + iPixel].y *= normal_fudge;
                pOutImageIQ[columnOffset + ir].x /= normCoef;
		        pOutImageIQ[columnOffset + ir].y /= normCoef;
		        pOutImageIQ[columnOffset + ir].x *= normal_fudge;
		        pOutImageIQ[columnOffset + ir].y *= normal_fudge;
		    }
		}
		else
		{
		    #pragma unroll 2
		    for (int iPixel = 0; iPixel < g_bf_params.m_nbLinesPerRecon; iPixel++)
		    {
                int columnOffset = (ix*g_bf_params.m_nbLinesPerRecon+iPixel) * g_bf_params.m_nbPixelsPerLine;
// 		        pOutImageIQ[columnOffset + ir + iPixel].x /= sqrt(normCoef);
// 		        pOutImageIQ[columnOffset + ir + iPixel].y /= sqrt(normCoef);
// 		        pOutImageIQ[columnOffset + ir + iPixel].x *= normal_fudge;
// 		        pOutImageIQ[columnOffset + ir + iPixel].y *= normal_fudge;
                pOutImageIQ[columnOffset + ir].x /= sqrt(normCoef);
		        pOutImageIQ[columnOffset + ir].y /= sqrt(normCoef);
		        pOutImageIQ[columnOffset + ir].x *= normal_fudge;
		        pOutImageIQ[columnOffset + ir].y *= normal_fudge;
		    }
		}
	    }
    }
}