///////////////////////////////////////////////////////
//                                                   //
//  Copyright © 2005-2011 by SuperSonic Imagine, SA  //
//         Confidential - All Right Reserved         //
//                                                   //
///////////////////////////////////////////////////////

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

const int NB_MAX_ANGLES = 81;
const int NB_ACQ_CHANNELS = 128 ;
const int HALF_MAX_APERTURE = 64 ;
const int HALF_MIN_APERTURE = 2;

struct SBeamformParameters
{
    float m_speedOfSound;
    float m_ProbeRadius;
    float m_invSpeedOfSound;
    float m_piezoPitch;
    float m_demodFreq;
    float m_rOrigin;
    float m_fNumber;
    float m_peakDelay;
    float m_xSources[NB_MAX_ANGLES];
    float m_zSources[NB_MAX_ANGLES];
//    float m_receiveAngle;
    float m_pixelPitch;
    float m_linePitch;
    float m_lambda;
	float m_xOrigin;
    float m_timeOrigin[NB_MAX_ANGLES];
    int   m_nbPiezos;
    int   m_nbAngles;
    int   m_channelOffset;
    int   m_firstSample;
    int   m_nbSamples;
    int   m_nbPixelsPerLine;
    int   m_nbRecon;
    int   m_firstAcquiredChannel;
    int   m_normMode;
    int   m_synthAcq;
    int   m_frame_per_frame;
    int   m_usegpu;
    int   m_nbImages;
};

void cudaIQBeamform(
    const SBeamformParameters * inParams,
    const short * inRF,
    const int bufferSize,
    float2 * outImageIQ);

void iq_beamform_kernel_cpu(const short* pInRF,const int bufferSize, float2* pOutImageIQ, SBeamformParameters g_bf_params);
