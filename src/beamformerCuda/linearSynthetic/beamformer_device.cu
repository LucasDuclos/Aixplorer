///////////////////////////////////////////////////////
//                                                   //
//  Copyright © 2005-2011 by SuperSonic Imagine, SA  //
//         Confidential - All Right Reserved         //
//                                                   //
///////////////////////////////////////////////////////
//
// libssip / src/ cuda / beamformer_device.cu
//

__device__ __constant__ SBeamformParameters g_bf_params;

__device__ float sinTable(int iAngle)
{
    return sinf(g_bf_params.m_angles[iAngle]) * g_bf_params.m_invSpeedOfSound;
}

__device__ float cosTable(int iAngle)
{
    return cosf(g_bf_params.m_angles[iAngle]) * g_bf_params.m_invSpeedOfSound;
}

__device__ float tanTable(int iAngle)
{
    return tanf(g_bf_params.m_angles[iAngle]);
}

__device__ float apod(int n, int apertureLength)
{
    return (0.54f - 0.46f * cosf(2.0f * (float)(M_PI) * (float)(n+1) / ((float)apertureLength - 1.0f)));
}