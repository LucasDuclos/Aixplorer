///////////////////////////////////////////////////////
//                                                   //
//  Copyright © 2005-2008 by SuperSonic Imagine, SA  //
//         Confidential - All Right Reserved         //
//                                                   //
///////////////////////////////////////////////////////
#pragma once

/// CompLut
/// @param Probe : double[3] to describe probe
/// @param Region, double[10] to describe imaged area
/// @param ReconInfo : double[15] to describe the way to reconstruct
/// @param ReconMux : double
/// @param ReconAbsc : double to describe abscissa
/// @param ReconDelay : double
/// @param Lut : short[LutSize] to contain computed Lut
/// @return 0 if ok; !=0 means problem
int CompLut(
        double *Probe, double *Region, double *ReconInfo,
        int *ReconMux, int *ReconAbsc, double * ReconDelay,
        short *Lut);
