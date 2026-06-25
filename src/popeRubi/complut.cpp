
#include "CompLutSplitted.h"

int CompLut(double *Probe, double *AcqInfo, double *ReconInfo, int *ReconMux, int *ReconAbsc, double * ReconDelay, short *Lut)
{
    if(Probe[3] == 0.0)
    {
        if (ReconInfo[1]==10 || ReconInfo[1] == 11)
        {
            return CompLutLinearPhotoacoustics(Probe,AcqInfo,ReconInfo,ReconMux,ReconAbsc,ReconDelay,Lut);
        }
        else
        {
            return CompLutLinear(Probe,AcqInfo,ReconInfo,ReconMux,ReconAbsc,ReconDelay,Lut);
        }
    }
    if(Probe[3] < 0.0)
    {
    	return CompLutPhased(Probe,AcqInfo,ReconInfo,ReconMux,ReconAbsc,ReconDelay,Lut);

    }
    else
    {
        if (ReconInfo[1] == 0)
    	{
            return CompLutCurvedSteering(Probe,AcqInfo,ReconInfo,ReconMux,ReconAbsc,ReconDelay,Lut);
    	}
		else
    	{
            return CompLut_FlatSynthetic_Curved(Probe,AcqInfo,ReconInfo,ReconMux,ReconAbsc,ReconDelay,Lut);
    	}
    }
}


