///////////////////////////////////////////////////////
//                                                   //
//  Copyright © 2005-2011 by SuperSonic Imagine, SA  //
//         Confidential - All Right Reserved         //
//                                                   //
///////////////////////////////////////////////////////

#include <stdio.h>
#include <cuda.h>


////////////////////////////////////////////////////////////////////////////////
// Useful CUDA error checking macros (adapted from similar macros in CUDA SDK cutil.h)
////////////////////////////////////////////////////////////////////////////////
#define DEBUG_CUDA
#ifdef DEBUG_CUDA

#define CUDA_SAFE_CALL_NO_SYNC( call ) do {                                  \
    cudaError err = call;                                                    \
    if (cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString(err));                \
    } } while (0)

#define CUDA_SAFE_CALL( call ) do {                                          \
    CUDA_SAFE_CALL_NO_SYNC(call);                                            \
    cudaError err = cudaThreadSynchronize();                                 \
    if (cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString(err));                \
    } } while (0)

#define CUDA_CHECK_ERROR( errorMessage ) do {                                \
    cudaError_t err = cudaGetLastError();                                    \
    if (cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString(err));  \
    }                                                                        \
    err = cudaThreadSynchronize();                                           \
    if (cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString(err));  \
    } } while (0)

#define CUFFT_SAFE_CALL( call ) do {                                         \
    cufftResult err = call;                                                  \
    if( CUFFT_SUCCESS != err) {                                              \
        fprintf(stderr, "CUFFT error in file '%s' in line %i.\n",            \
                __FILE__, __LINE__);                                         \
    } } while (0)


#else

#define CUDA_SAFE_CALL_NO_SYNC( call ) call
#define CUDA_SAFE_CALL( call ) call
#define CUDA_CHECK_ERROR( errorMessage )
#define CUFFT_SAFE_CALL( call ) call

#endif // DEBUG

// Round a / b to nearest higher integer value
#define IDIVUP(a, b) ((a % b != 0) ? (a / b + 1) : (a / b))

// Fast (for 1.x Cuda device) integer multiplication
#if __CUDA_ARCH__ >= 200
#define IMUL(a, b) ((a)*(b))
#else
#define IMUL(a, b) __mul24(a, b)
#endif