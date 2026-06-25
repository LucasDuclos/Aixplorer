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
        return a ;
    else
        return b ;
}

inline int max(int a, int b)
{
    if(a>b)
        return a ;
    else
        return b ;
}


float apod_c(int n, int apertureLength)
{
    return (0.54f - 0.46f * cos(2.0f * (float)(M_PI) * (float)(n+1) / ((float)apertureLength - 1.0f)));
}

float sinTable(int iAngle,SBeamformParameters g_bf_params)
{
#ifdef USE_LUTS
    return g_bf_sinTable[iAngle];
#else
    return sin(g_bf_params.m_angles[iAngle]) * g_bf_params.m_invSpeedOfSound;
#endif
}

float cosTable(int iAngle,SBeamformParameters g_bf_params)
{
#ifdef USE_LUTS
    return g_bf_cosTable[iAngle];
#else
    return cos(g_bf_params.m_angles[iAngle]) * g_bf_params.m_invSpeedOfSound;
#endif
}

float tanTable(int iAngle,SBeamformParameters g_bf_params)
{
#ifdef USE_LUTS
    return g_bf_sinTable[iAngle]/g_bf_cosTable[iAngle];
#else
    return tan(g_bf_params.m_angles[iAngle]);
#endif
}

void iq_beamform_kernel_cpu(const short* pInRF,const int bufferSize, float2* pOutImageIQ, SBeamformParameters g_bf_params)
{
    // output pixel coords
    
    for(int ix=0 ; ix< g_bf_params.m_nbRecon ; ix++)
    {
        for(int iz=0 ; iz<g_bf_params.m_nbPixelsPerLine ; iz++)
        {

		float Z = (float)iz * g_bf_params.m_pixelPitch * cos(g_bf_params.m_receiveAngle) + g_bf_params.m_zOrigin;
		float aperture = Z/g_bf_params.m_fNumber;
		int halfAperture = (int) round(0.5f * aperture / g_bf_params.m_piezoPitch);
		halfAperture = min(halfAperture, HALF_MAX_APERTURE);
		halfAperture = max(halfAperture, HALF_MIN_APERTURE);

		//int outputOffset = nbLines * ix + iz * g_bf_params.m_nbLinesPerRecon;
		int columnOffset = ix * g_bf_params.m_nbLinesPerRecon * g_bf_params.m_nbPixelsPerLine;

		float refDelay = 2.0f * Z * g_bf_params.m_invSpeedOfSound + g_bf_params.m_peakDelay;
		
		const float w0 = -1.0f/6.0f;
		const float w1 = 0.5f;
		const float w2 = -0.5f;
		const float w3 = 1.0f/6.0f;
		const float eps = 0.000001f;
		const float normal_fudge = 1.0f;

		int piezoRef = min(g_bf_params.m_firstRecon + (int)round((float)ix*g_bf_params.m_reconPitch/g_bf_params.m_piezoPitch), g_bf_params.m_nbPiezos);
		int firstChannel = max(g_bf_params.m_firstAcquiredChannel, piezoRef - halfAperture);
		int lastChannel  = min(g_bf_params.m_firstAcquiredChannel + NB_ACQ_CHANNELS * (1 + g_bf_params.m_synthAcq), piezoRef + halfAperture);
		float xPixelLeft = ((float)g_bf_params.m_firstRecon-(float)(g_bf_params.m_nbPiezos/2) + g_bf_params.m_xOrigin)*g_bf_params.m_piezoPitch+g_bf_params.m_reconPitch*ix
		                 - g_bf_params.m_reconPitch/2.0f + Z * tan(g_bf_params.m_receiveAngle);
		
        
        int firstAngle = 0;
		int lastAngle = g_bf_params.m_nbAngles;
		
		float normCoef = eps;
		float IQ[2];
		int offset;

		for (int channel = firstChannel; channel < lastChannel; channel++)
		{
	#ifdef USE_CHANNEL_OFFSET
		    offset = g_bf_chanOffsets[channel];
	#else
		    offset = (channel % NB_ACQ_CHANNELS) * g_bf_params.m_channelOffset + g_bf_params.m_synthAcq * (channel / NB_ACQ_CHANNELS) * g_bf_params.m_nbSamples;
		   
	#endif
            
		    float apodCoef = apod_c(channel - firstChannel, lastChannel - firstChannel+1);
		    float xPiezo = (float)(channel - g_bf_params.m_nbPiezos/2) * g_bf_params.m_piezoPitch;
		    
		    for (int iAngle = firstAngle; iAngle < lastAngle; iAngle++)
		    {
		        const short* channelStartPtr = pInRF + offset + (1+g_bf_params.m_synthAcq)*iAngle*g_bf_params.m_nbSamples;
		        float angleApod = apod_c(iAngle - firstAngle, lastAngle - firstAngle+1);

		        #pragma unroll 2
		        for (int iPixel = 0; iPixel < g_bf_params.m_nbLinesPerRecon; iPixel++)
		        {
		            float X = xPixelLeft + iPixel * g_bf_params.m_linePitch;
		            float forwardDelay = X * sinTable(iAngle,g_bf_params) + Z * cosTable(iAngle,g_bf_params) + g_bf_params.m_timeOrigin[iAngle] + g_bf_params.m_peakDelay;
		            float returnDelay  = sqrt(Z * Z + (X - xPiezo) * (X - xPiezo)) * g_bf_params.m_invSpeedOfSound;
		            float timeDelay = forwardDelay + returnDelay;
                   
		            int delay = (int) floor(timeDelay * g_bf_params.m_demodFreq) - g_bf_params.m_firstSample/2;
                    
		            float deltaDelay = 2.0f * (float)(M_PI) * (timeDelay - refDelay) * g_bf_params.m_demodFreq;
		            float phi = timeDelay * g_bf_params.m_demodFreq - floor(timeDelay * g_bf_params.m_demodFreq);
	
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
		            
                    pOutImageIQ[columnOffset + iz + iPixel].x += (IQ[0]*cos(deltaDelay) + IQ[1]*sin(deltaDelay))*apodCoef;
		            pOutImageIQ[columnOffset + iz + iPixel].y += (-IQ[0]*sin(deltaDelay) + IQ[1]*cos(deltaDelay))*apodCoef;
		        }
		    }
		    normCoef += apodCoef * apodCoef;
		}

		if (g_bf_params.m_normMode < 100)
		{
		    #pragma unroll 2
		    for (int iPixel = 0; iPixel < g_bf_params.m_nbLinesPerRecon; iPixel++)
		    {
		        pOutImageIQ[columnOffset + iz + iPixel].x /= normCoef;
		        pOutImageIQ[columnOffset + iz + iPixel].y /= normCoef;
		        pOutImageIQ[columnOffset + iz + iPixel].x *= normal_fudge;
		        pOutImageIQ[columnOffset + iz + iPixel].y *= normal_fudge;
		    }
		}
		else
		{
		    #pragma unroll 2
		    for (int iPixel = 0; iPixel < g_bf_params.m_nbLinesPerRecon; iPixel++)
		    {
		        pOutImageIQ[columnOffset + iz + iPixel].x /= sqrt(normCoef);
		        pOutImageIQ[columnOffset + iz + iPixel].y /= sqrt(normCoef);
		        pOutImageIQ[columnOffset + iz + iPixel].x *= normal_fudge;
		        pOutImageIQ[columnOffset + iz + iPixel].y *= normal_fudge;
		    }
		}
	    }
    }
}