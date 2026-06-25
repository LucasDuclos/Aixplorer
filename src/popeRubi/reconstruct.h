///////////////////////////////////////////////////////
//                                                   //
//  Copyright © 2005-2008 by SuperSonic Imagine, SA  //
//         Confidential - All Right Reserved         //
//                                                   //
///////////////////////////////////////////////////////
//
// sp / reconstruct.h
//
#pragma once

#ifdef __linux__
    #include "pthread.h"
#elif defined _WIN32 || defined _WIN64
    #include <windows.h>
    #include <stdlib.h>
    #include <string.h>
    #include <stdio.h>
    #include <conio.h>
    #include <process.h>
#else
    #error "unknown platform"
#endif

enum EReconMethod
{
    E_INTRINSICS =0,
    E_ASM
};

enum EDataFormat
{
    E_RUBI =0,
    E_MATLAB
};


struct ArgReconstruct
{
short * Lut ;
short * RcvDataShort ;
float * IQPix ;
float * IQPixImag ;
float * IQbuffer ;
int     dataFormat ;
};

namespace sp
{
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////    
    void * reconstruct(void  *arg);

}
