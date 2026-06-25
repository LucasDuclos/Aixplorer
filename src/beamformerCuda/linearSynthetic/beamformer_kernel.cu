///////////////////////////////////////////////////////
//                                                   //
//  Copyright © 2005-2011 by SuperSonic Imagine, SA  //
//         Confidential - All Right Reserved         //
//                                                   //
///////////////////////////////////////////////////////
//
// cuda/beamformer_kernel.cu
//

#include "beamformerCu.h"

__global__ void iq_beamform_kernel(const short* pInRF,const int bufferSize, float2* pOutImageIQ)
{
    
    // output pixel coords
    int ix = IMUL(blockDim.x, blockIdx.x) + threadIdx.x;  // recon
    int iy = IMUL(blockDim.y, blockIdx.y) + threadIdx.y;  // sample (depth)
    if(ix<g_bf_params.m_nbRecon && iy < g_bf_params.m_nbPixelsPerLine)
    {

		float Z = (float)iy * g_bf_params.m_pixelPitch * cos(g_bf_params.m_receiveAngle) + g_bf_params.m_zOrigin;
		float aperture = Z/g_bf_params.m_fNumber;
		int halfAperture = (int) roundf(0.5f * aperture / g_bf_params.m_piezoPitch);
		halfAperture = min(halfAperture, HALF_MAX_APERTURE);
		halfAperture = max(halfAperture, HALF_MIN_APERTURE);
		
		float refDelay = 2.0f * Z * g_bf_params.m_invSpeedOfSound + g_bf_params.m_peakDelay;
		
		const float w0 = -1.0f/6.0f;
		const float w1 = 0.5f;
		const float w2 = -0.5f;
		const float w3 = 1.0f/6.0f;
		const float eps = 0.000001f;
		const float normal_fudge = 1.0f;

		int piezoRef = min(g_bf_params.m_firstRecon + (int)roundf((float)ix*g_bf_params.m_reconPitch/g_bf_params.m_piezoPitch), g_bf_params.m_nbPiezos);
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

		    offset = (channel % NB_ACQ_CHANNELS) * g_bf_params.m_channelOffset + g_bf_params.m_synthAcq * (channel / NB_ACQ_CHANNELS) * g_bf_params.m_nbSamples;
		   
		    float apodCoef = apod(channel - firstChannel, lastChannel - firstChannel+1);
		    float xPiezo = (float)(channel - g_bf_params.m_nbPiezos/2) * g_bf_params.m_piezoPitch;
		    
		    for (int iAngle = firstAngle; iAngle < lastAngle; iAngle++)
		    {
		        const short* channelStartPtr = pInRF + offset + (1+g_bf_params.m_synthAcq)*iAngle*g_bf_params.m_nbSamples;
		        float angleApod = apod(iAngle - firstAngle, lastAngle - firstAngle+1);

		        for (int iPixel = 0; iPixel < g_bf_params.m_nbLinesPerRecon; iPixel++)
		        {
		            float X = xPixelLeft + iPixel * g_bf_params.m_linePitch;
		            float forwardDelay = X * sinTable(iAngle) + Z * cosTable(iAngle) + g_bf_params.m_timeOrigin[iAngle] + g_bf_params.m_peakDelay;
		            float returnDelay  = sqrtf(Z * Z + (X - xPiezo) * (X - xPiezo)) * g_bf_params.m_invSpeedOfSound;
		            float timeDelay = forwardDelay + returnDelay;
                   
		            int delay = (int) floorf(timeDelay * g_bf_params.m_demodFreq) - g_bf_params.m_firstSample/2;

                if(delay < 1)
                {
                    delay = 1 ;
                }                    
                
                if( 2*delay-2+1 + offset >= bufferSize )
                {
                    delay = 1 ;
                }
                    
		            float deltaDelay = 2.0f * (float)(M_PI) * (timeDelay - refDelay) * g_bf_params.m_demodFreq;
		            float phi = timeDelay * g_bf_params.m_demodFreq - floorf(timeDelay * g_bf_params.m_demodFreq);
	
		            float X0 = phi+1.0f;
		            float X1 = phi+eps;
		            float X2 = phi-1.0f;
		            float X3 = phi-2.0f;
		            
		            // expand some terms and eliminate divides
		            float l0 = X1*X2*X3*w0;
		            float l1 = X0*X2*X3*w1;
		            float l2 = X0*X1*X3*w2;
		            float l3 = X0*X1*X2*w3;
		            IQ[0] = l0*(float)channelStartPtr[2*delay-2] + l1*(float)channelStartPtr[2*delay] + l2*(float)channelStartPtr[2*delay+2] + l3*(float)channelStartPtr[2*delay+4];
		            IQ[1] = l0*(float)channelStartPtr[2*delay-2+1] + l1*(float)channelStartPtr[2*delay+1] + l2*(float)channelStartPtr[2*delay+2+1] + l3*(float)channelStartPtr[2*delay+4+1];

                    int columnOffset = (ix*g_bf_params.m_nbLinesPerRecon + iPixel)  * g_bf_params.m_nbPixelsPerLine;
		            pOutImageIQ[columnOffset + iy].x += (IQ[0]*cos(deltaDelay) + IQ[1]*sin(deltaDelay))*apodCoef;
		            pOutImageIQ[columnOffset + iy].y += (-IQ[0]*sin(deltaDelay) + IQ[1]*cos(deltaDelay))*apodCoef;
		        }
		    }
		    normCoef += apodCoef * apodCoef;
		}

		if (g_bf_params.m_normMode < 100)
		{
		    #pragma unroll 2
		    for (int iPixel = 0; iPixel < g_bf_params.m_nbLinesPerRecon; iPixel++)
		    {
                int columnOffset = (ix*g_bf_params.m_nbLinesPerRecon + iPixel)  * g_bf_params.m_nbPixelsPerLine;
		        pOutImageIQ[columnOffset + iy].x /= normCoef;
		        pOutImageIQ[columnOffset + iy].y /= normCoef;
		        pOutImageIQ[columnOffset + iy].x *= normal_fudge;
		        pOutImageIQ[columnOffset + iy].y *= normal_fudge;
		    }
		}
		else
		{
		    #pragma unroll 2
		    for (int iPixel = 0; iPixel < g_bf_params.m_nbLinesPerRecon; iPixel++)
		    {
                int columnOffset = (ix*g_bf_params.m_nbLinesPerRecon + iPixel)  * g_bf_params.m_nbPixelsPerLine;
		        pOutImageIQ[columnOffset + iy].x /= sqrtf(normCoef);
		        pOutImageIQ[columnOffset + iy].y /= sqrtf(normCoef);
		        pOutImageIQ[columnOffset + iy].x *= normal_fudge;
		        pOutImageIQ[columnOffset + iy].y *= normal_fudge;
		    }
		}
    }
}
