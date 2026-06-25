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

        float thetaFiring = g_bf_params.m_firstFiringAngle + ix*g_bf_params.m_firingAnglePitch;
        float RfromApex = (float)iy * g_bf_params.m_pixelPitch + g_bf_params.m_rOrigin;
        float aperture = RfromApex/g_bf_params.m_fNumber;
		int halfAperture = (int) roundf(0.5f * aperture / g_bf_params.m_piezoPitch);
        int firstChannel0 = g_bf_params.m_nbPiezos/2 - halfAperture; 
        int lastChannel0 = g_bf_params.m_nbPiezos/2 + halfAperture;
		
		const float w0 = -1.0f/6.0f;
		const float w1 = 0.5f;
		const float w2 = -0.5f;
		const float w3 = 1.0f/6.0f;
		const float eps = 0.000001f;
		const float normal_fudge = 1.0f;

		int firstChannel = 0;
		int lastChannel  = g_bf_params.m_nbPiezos;
				
		float normCoef = eps;
		float IQ[2];
		
		for (int channel = firstChannel; channel < lastChannel; channel++)
		{
	        int offset = (channel % NB_ACQ_CHANNELS) * g_bf_params.m_channelOffset + ix * g_bf_params.m_nbSamples;
		   
		    float apodCoef = apod(channel - firstChannel0, lastChannel0 - firstChannel0+1);
		    float xPiezo = (float)(channel +0.5f - g_bf_params.m_nbPiezos/2) * g_bf_params.m_piezoPitch;
            float zPiezo = g_bf_params.m_zApex;
		    		    
            const short* channelStartPtr = pInRF + offset;
            
            for (int iPixel = 0; iPixel < g_bf_params.m_nbLinesPerRecon; iPixel++)
            {
                int columnOffset = (ix*g_bf_params.m_nbLinesPerRecon+iPixel) * g_bf_params.m_nbPixelsPerLine;
                float thetaLine = thetaFiring + iPixel * g_bf_params.m_linePitch + g_bf_params.m_thetaOrigin;
                float costhetaLine = cosf(thetaLine);
                float sinthetaLine = sinf(thetaLine);
                float RfromSkinLine = RfromApex - zPiezo/costhetaLine;
                float X = RfromApex*sinthetaLine;
                float Z = RfromApex*costhetaLine;
                
                float forwardDelay = RfromSkinLine*g_bf_params.m_invSpeedOfSound;
                float returnDelay  = sqrtf((Z-zPiezo) * (Z-zPiezo) + (X - xPiezo) * (X - xPiezo)) * g_bf_params.m_invSpeedOfSound;
                float timeDelay = forwardDelay + returnDelay + g_bf_params.m_timeOrigin[ix];

                int delay = (int) floorf(timeDelay * g_bf_params.m_demodFreq) - g_bf_params.m_firstSample/2;
                
                if(delay < 0)
                {
                    delay = 1 ;
                }                    
                
                if( 2*delay-2+1 + offset >= bufferSize )
                {
                    delay = 1 ;
                }
                float deltaDelay = 2.0f * (float)(M_PI) * (timeDelay) * g_bf_params.m_demodFreq;
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
                IQ[0] = l0*channelStartPtr[2*delay-2] + l1*channelStartPtr[2*delay] + l2*channelStartPtr[2*delay+2] + l3*channelStartPtr[2*delay+4];
                IQ[1] = l0*channelStartPtr[2*delay-2+1] + l1*channelStartPtr[2*delay+1] + l2*channelStartPtr[2*delay+2+1] + l3*channelStartPtr[2*delay+4+1];

                pOutImageIQ[columnOffset + iy].x += (IQ[0]*cosf(deltaDelay) + IQ[1]*sinf(deltaDelay))*apodCoef;
                pOutImageIQ[columnOffset + iy].y += (-IQ[0]*sinf(deltaDelay) + IQ[1]*cosf(deltaDelay))*apodCoef;
            }

		    normCoef += apodCoef * apodCoef;
		}
        
		if (g_bf_params.m_normMode < 100)
		{
#pragma unroll 2
		    for (int iPixel = 0; iPixel < g_bf_params.m_nbLinesPerRecon; iPixel++)
		    {
                int columnOffset = (ix*g_bf_params.m_nbLinesPerRecon+iPixel) * g_bf_params.m_nbPixelsPerLine;
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
                int columnOffset = (ix*g_bf_params.m_nbLinesPerRecon+iPixel) * g_bf_params.m_nbPixelsPerLine;
                pOutImageIQ[columnOffset + iy].x /= sqrtf(normCoef);
		        pOutImageIQ[columnOffset + iy].y /= sqrtf(normCoef);
		        pOutImageIQ[columnOffset + iy].x *= normal_fudge;
		        pOutImageIQ[columnOffset + iy].y *= normal_fudge;
		    }
		}
    }
}
