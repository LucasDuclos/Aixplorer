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
        const int imgOffset = ix* g_bf_params.m_nbPixelsPerLine + iy;

		float Ro = (float)iy * g_bf_params.m_pixelPitch + g_bf_params.m_rOrigin;
		float aperture = Ro/g_bf_params.m_fNumber;
		int halfAperture = (int) roundf(0.5f * aperture / g_bf_params.m_piezoPitch);
		halfAperture = min(halfAperture, HALF_MAX_APERTURE);
		halfAperture = max(halfAperture, HALF_MIN_APERTURE);
		
		float refDelay = 0;
		
		const float w0 = -1.0f/6.0f;
		const float w1 = 0.5f;
		const float w2 = -0.5f;
		const float w3 = 1.0f/6.0f;
		const float eps = 0.000001f;
		const float normal_fudge = 1.0f;

		int piezoRef = min((int)(g_bf_params.m_xOrigin + 0.5f + ( ix*g_bf_params.m_linePitch ) /g_bf_params.m_piezoPitch ), g_bf_params.m_nbPiezos);
		int firstChannel = max(g_bf_params.m_firstAcquiredChannel, piezoRef - halfAperture);
		int lastChannel  = min(g_bf_params.m_firstAcquiredChannel + NB_ACQ_CHANNELS * (1 + g_bf_params.m_synthAcq), piezoRef + halfAperture);
		
        float xPixel = ((float)g_bf_params.m_xOrigin-(float)(g_bf_params.m_nbPiezos/2))*g_bf_params.m_piezoPitch + ix*g_bf_params.m_linePitch;
		
        lastChannel = min(lastChannel,g_bf_params.m_nbPiezos);

        int firstAngle = 0;
		int lastAngle = g_bf_params.m_nbAngles;

        if(g_bf_params.m_frame_per_frame>=0 && g_bf_params.m_frame_per_frame<lastAngle)
        {
            firstAngle = g_bf_params.m_frame_per_frame;
            lastAngle = g_bf_params.m_frame_per_frame + 1;
        }
		
		float normCoef = eps;
		float IQ[2];
		int offset;

		for (int channel = firstChannel; channel < lastChannel; channel++)
		{

		    offset = (channel % NB_ACQ_CHANNELS) * g_bf_params.m_channelOffset + g_bf_params.m_synthAcq * (channel / NB_ACQ_CHANNELS) * g_bf_params.m_nbSamples;
		   
		    float apodCoef = apod(channel - firstChannel, lastChannel - firstChannel+1);
		    
            float xPiezo = g_bf_params.m_ProbeRadius*sinf((float)(channel - g_bf_params.m_nbPiezos/2) * g_bf_params.m_piezoPitch/g_bf_params.m_ProbeRadius);
            float zPiezo = g_bf_params.m_ProbeRadius*cosf((float)(channel - g_bf_params.m_nbPiezos/2) * g_bf_params.m_piezoPitch/g_bf_params.m_ProbeRadius);

		    for (int iAngle = firstAngle; iAngle < lastAngle; iAngle++)
		    {
		        const short* channelStartPtr = pInRF + offset + (1+g_bf_params.m_synthAcq)*iAngle*g_bf_params.m_nbSamples;

                float aPolar = xPixel/g_bf_params.m_ProbeRadius;
                float rPolar = Ro + g_bf_params.m_ProbeRadius;
                float X = rPolar*sinf(aPolar);
                float Z = rPolar*cosf(aPolar);
                
                float xSource = g_bf_params.m_xSources[iAngle];
                float zSource = g_bf_params.m_zSources[iAngle];

                float forwardDelay = sqrtf( (X - xSource) * (X - xSource) + (Z - zSource) * (Z - zSource) ) * g_bf_params.m_invSpeedOfSound;
                float returnDelay  = sqrtf( (Z - zPiezo)  * (Z - zPiezo)  + (X - xPiezo)  * (X - xPiezo)  ) * g_bf_params.m_invSpeedOfSound;
                
                float timeDelay = forwardDelay + returnDelay - g_bf_params.m_timeOrigin[iAngle];

                int delay = (int) floorf(timeDelay * g_bf_params.m_demodFreq) - g_bf_params.m_firstSample/2;

                if(delay < 1)
                {
                    delay = 1 ;
                }                    

                if( 2*delay-2+1 + offset >= bufferSize )
                {
                    delay = 1 ;
                }
                float deltaDelay = 2.0f * (float)(M_PI) * (timeDelay -refDelay) * g_bf_params.m_demodFreq;
                float phi = timeDelay * g_bf_params.m_demodFreq - floorf(timeDelay * g_bf_params.m_demodFreq);

                float X0 = phi+1.0f;
                float X1 = phi+eps;
                float X2 = phi-1.0f;
                float X3 = phi-2.0f;

                float X0X1 = X0*X1 ;
                float X2X3 = X2*X3 ;

                // expand some terms and eliminate divides
                float l0 = X1*X2X3*w0;
                float l1 = X0*X2X3*w1;
                float l2 = X0X1*X3*w2;
                float l3 = X0X1*X2*w3;

                IQ[0] = l0*(float)channelStartPtr[2*delay-2] + l1*(float)channelStartPtr[2*delay] + l2*(float)channelStartPtr[2*delay+2] + l3*(float)channelStartPtr[2*delay+4];
                IQ[1] = l0*(float)channelStartPtr[2*delay-2+1] + l1*(float)channelStartPtr[2*delay+1] + l2*(float)channelStartPtr[2*delay+2+1] + l3*(float)channelStartPtr[2*delay+4+1];
                
                float cos1 = cosf(deltaDelay);
                float sin1 = sinf(deltaDelay);
                
                if(iAngle==0 && channel == firstChannel)
                {
                    pOutImageIQ[imgOffset].x = (IQ[0]*cos1 + IQ[1]*sin1)*apodCoef;
                    pOutImageIQ[imgOffset].y = (-IQ[0]*sin1 + IQ[1]*cos1)*apodCoef;
                }
                else
                {
                    pOutImageIQ[imgOffset].x += (IQ[0]*cos1 + IQ[1]*sin1)*apodCoef;
                    pOutImageIQ[imgOffset].y += (-IQ[0]*sin1 + IQ[1]*cos1)*apodCoef;
                }
		    }
		    normCoef += apodCoef * apodCoef;
		}

		if (g_bf_params.m_normMode < 100)
		{
            pOutImageIQ[imgOffset].x /= normCoef;
            pOutImageIQ[imgOffset].y /= normCoef;
            pOutImageIQ[imgOffset].x *= normal_fudge;
            pOutImageIQ[imgOffset].y *= normal_fudge;
    }
		else
		{
            pOutImageIQ[imgOffset].x /= sqrtf(normCoef);
            pOutImageIQ[imgOffset].y /= sqrtf(normCoef);
            pOutImageIQ[imgOffset].x *= normal_fudge;
            pOutImageIQ[imgOffset].y *= normal_fudge;
        }
    }
}
