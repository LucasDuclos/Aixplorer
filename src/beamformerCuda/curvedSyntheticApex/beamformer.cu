///////////////////////////////////////////////////////
//                                                   //
//  Copyright © 2005-2011 by SuperSonic Imagine, SA  //
//         Confidential - All Right Reserved         //
//                                                   //
///////////////////////////////////////////////////////
//
// libssip / src/ cuda / beamformer.cu
//

#include "beamformer.h"
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>

// The IQ beamform works best with a smaller block size (8,8) rather than (16,16).
// TODO: Investigate this.
#define NB_THREADS_X 8
#define NB_THREADS_Y 8

// Uncomment to use lookup tables for channel offset to eliminate the if-test (and divergent branches) in the Cuda kernel
//#define USE_CHANNEL_OFFSET

// Uncomment to use eps check in V2 kernel
//#define USE_BRANCH


// Get common device functions and data.
#include "beamformer_device.cu"
#include "beamformer_kernel.cu"

void cudaIQBeamform(
    const SBeamformParameters * inParams,
    const short * inRF,
    const int bufferSize,
    float2 * outImageIQ)
{
    short * pRFgpu ; // pointer to the Rfdata in the gpu
    float2 * pIQOutDatagpu ; // pointer to the output data in the gpu
    
    if(inParams->m_usegpu)
    {
        cudaDeviceReset();
    }

    // Increase the L1 cache size to 48kB (at the expense of global shared memory).
    cudaFuncSetCacheConfig(iq_beamform_kernel, cudaFuncCachePreferL1);


    const int imgSize = inParams->m_nbPixelsPerLine * inParams->m_nbRecon ;
   
    // allocation of the gpu ram

    CUDA_SAFE_CALL( cudaMalloc( (void**)&pRFgpu, bufferSize*sizeof(short) ) );
    CUDA_SAFE_CALL( cudaMalloc( (void**)&pIQOutDatagpu, imgSize*sizeof(float2) ) );
    CUDA_SAFE_CALL( cudaMemset( pIQOutDatagpu, 0 , imgSize*sizeof(float2) ) );

    // copy bf_params from host to device
    CUDA_SAFE_CALL( cudaMemcpyToSymbol(g_bf_params, inParams, sizeof(SBeamformParameters)) );

    // copy RF data from host to device
    CUDA_SAFE_CALL( cudaMemcpy( pRFgpu, inRF , bufferSize*sizeof(short), cudaMemcpyHostToDevice));

    dim3 block( NB_THREADS_X, NB_THREADS_Y );
    dim3 grid( IDIVUP(inParams->m_nbRecon, block.x), IDIVUP(inParams->m_nbPixelsPerLine, block.y) );
    
    for(int iImage = 0 ; iImage < inParams->m_nbImages ; iImage++)
    {
    
        iq_beamform_kernel<<<grid, block>>>(&pRFgpu[iImage * inParams->m_nbSamples * (1 + inParams->m_synthAcq) * inParams->m_nbAngles], bufferSize, pIQOutDatagpu);
        CUDA_CHECK_ERROR("iq_beamform_kernel");

        // copy final image from device to host
        CUDA_SAFE_CALL( cudaMemcpy( outImageIQ + iImage * imgSize, pIQOutDatagpu , imgSize*sizeof(float2), cudaMemcpyDeviceToHost));
    }

    // free GPU RAM
    CUDA_SAFE_CALL( cudaFree( pRFgpu ) );
    CUDA_SAFE_CALL( cudaFree( pIQOutDatagpu ) );
}