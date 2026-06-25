///////////////////////////////////////////////////////
//                                                   //
//  Copyright © 2005-2011 by SuperSonic Imagine, SA  //
//         Confidential - All Right Reserved         //
//                                                   //
///////////////////////////////////////////////////////

const int NB_MAX_SOURCES = 81;
const int NB_ACQ_CHANNELS = 128;
const int HALF_MAX_APERTURE = 64;
const int HALF_MIN_APERTURE = 2;

struct SBeamformParameters
{
    float m_speedOfSound;
    float m_invSpeedOfSound;
    float m_piezoPitch;
    float m_demodFreq;
    float m_rOrigin;
    float m_fNumber;
    float m_peakDelay;
    float m_xSource[NB_MAX_SOURCES];
    float m_zSource[NB_MAX_SOURCES];
    float m_zApex ;
    float m_linePitch;
    float m_pixelPitch;
    float m_lambda;
	float m_thetaOrigin;
    int   m_nbPiezos;
    int   m_nbSources;
    int   m_channelOffset;
    int   m_firstSample;
    int   m_nbSamples;
    int   m_nbLinesPerRecon;
    int   m_nbPixelsPerLine;
    int   m_nbRecon;
    int   m_normMode;
    int   m_synthAcq;
    int   m_frame_per_frame;
    int   m_idxTransmitToBeamform;
    int   m_usegpu;
};

void cudaIQBeamform(
    const SBeamformParameters * inParams,
    const short * inRF,
    const int bufferSize,
    float2 * outImageIQ, float & elapsedTimeMS);

void iq_beamform_kernel_cpu(const short* pInRF,const int bufferSize, float2* pOutImageIQ, SBeamformParameters g_bf_params);