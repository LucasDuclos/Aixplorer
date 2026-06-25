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

__device__ float apod(int n, int apertureLength)
{
    return (0.54f - 0.46f * cosf(2.0f * (float)(M_PI) * (float)(n+1) / ((float)apertureLength - 1.0f)));
}