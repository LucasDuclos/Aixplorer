#include "stdio.h"
#include "stdlib.h"
#include "mex.h"
#include "math.h"
#include "parameters.h"
#include "ParametersMethods.h"
#include "ScanConverterParameters.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))


void CompLut(int NR, int NTheta,float AngleOfFirstLine,float AngleStep,float FirstSampleAxialPosition,float RadialStep,float * OffSet,float * SizeOfImg,int Nx,int Nz,int NbVois,float SmoothDistance,float * XAxis,float * ZAxis,bool ComputeOffSetAndSize,
			float * LUT_coeff,short * LUT_X, short * LUT_Z ,short * LUT_Vois1, short * LUT_Vois2,
			int & LUTSIZE)
{
	float Sigma ;
	float dx,dz ;
	float Theta0 ;
	float * R0, *X0, *Z0 ;
	float Distance ;
	
	float SizeOfImageX, SizeOfImageZ ;
	float LUT_coeffSum ;

	int PositionLineUnit, PositionRadialSampleUnit ;
	int ClosiestLine, ClosiestAxialSample ;
	int SeekOffSet ;
	int k ;
	int indexVoisin ;

	float FinalSample_X, FinalSample_Z, FinalSample_R, FinalSample_theta ;
	float OffSetX, OffSetZ ;
	float temp1 , temp2 ;
	
	float * SINTHETA0, *COSTHETA0 ;
	float maxCos ;
	
	int NSample,NbVoisTot ;


	Sigma = 2*SmoothDistance ;
	Sigma *= Sigma ;

	NSample = NR*NTheta ;

	NbVoisTot = (2*NbVois) ; // Voisins
	NbVoisTot *=NbVoisTot ;

	X0 = (float*)malloc(NTheta*NR*sizeof(float));
	Z0 = (float*)malloc(NTheta*NR*sizeof(float));
	R0 = (float*)malloc(NTheta*NR*sizeof(float));
	
	SINTHETA0 = (float*)malloc(sizeof(float)*NTheta);
	COSTHETA0 = (float*)malloc(sizeof(float)*NTheta);
	
	// Trigonometric table
	maxCos = 0 ;
	for(int theta=0 ; theta<NTheta;theta++)
	{
		Theta0 = AngleOfFirstLine + theta*AngleStep ;
		SINTHETA0[theta] = sin(Theta0);
		COSTHETA0[theta] = cos(Theta0);
		if(COSTHETA0[theta]>maxCos)
			maxCos = COSTHETA0[theta] ;

	}

	// Offset and image size (computed if not passed)
	
	if(ComputeOffSetAndSize)
	{
		temp1 = (FirstSampleAxialPosition + RadialStep*(NR-1))*SINTHETA0[0] ;
		temp2 = FirstSampleAxialPosition*SINTHETA0[0] ;

		OffSetX = MIN(temp1,temp2);

		temp1 = FirstSampleAxialPosition+RadialStep*COSTHETA0[0] ;
		temp2 = FirstSampleAxialPosition*COSTHETA0[NTheta-1];

		OffSetZ = MIN(temp1,temp2);

		temp1 = (FirstSampleAxialPosition + RadialStep*(NR-1))*SINTHETA0[NTheta-1] ;
		temp2 = FirstSampleAxialPosition*SINTHETA0[NTheta-1] ;

		SizeOfImageX = MAX(temp1,temp2);
		SizeOfImageX -= OffSetX ;

		SizeOfImageZ = (FirstSampleAxialPosition + RadialStep*(NR-1))*maxCos ;

		SizeOfImageZ -= OffSetZ ;

		SizeOfImg[0] = SizeOfImageX ;
		SizeOfImg[1] = SizeOfImageZ ;
		OffSet[0] = OffSetX ;
		OffSet[1] = OffSetZ ;

	}
	else
	{
		SizeOfImageX = SizeOfImg[0] ;
		SizeOfImageZ = SizeOfImg[1] ;
		OffSetX = OffSet[0] ;
		OffSetZ = OffSet[1] ;

	}
	

	dx = SizeOfImageX/(Nx-1) ;
	dz = SizeOfImageZ/(Nz-1) ;

	//mexPrintf("dx = %g, dz = %g\n",dx,dz);

	// Cartesian coordinates table of input data

	for(int theta=0 ; theta<NTheta;theta++)
	{
		//Theta0 = AngleOfFirstLine + theta*AngleStep ;
		for(int r=0;r<NR;r++)
		{
			R0[theta+NTheta*r] = FirstSampleAxialPosition + r*RadialStep;
			X0[theta+NTheta*r] = R0[theta+NTheta*r]*SINTHETA0[theta];
			Z0[theta+NTheta*r] = R0[theta+NTheta*r]*COSTHETA0[theta];     			
		}
	}


	// Write Axis scale for display

	for(int x=0 ; x<Nx ; x++)
		XAxis[x] = x*dx + OffSetX ;
	
	for(int z=0 ; z<Nz ; z++)
	{
		ZAxis[z] = z*dz + OffSetZ ;
	}


	// first pass
	
	if(LUTSIZE==0)
		
		for(int x=0;x<Nx;x++)
		{
			FinalSample_X = XAxis[x] ;
			
			for(int z=0;z<Nz;z++)
			{
				FinalSample_Z = ZAxis[z] ;
				FinalSample_R = sqrt(FinalSample_X*FinalSample_X + FinalSample_Z*FinalSample_Z);
				FinalSample_theta = atan2(FinalSample_X,FinalSample_Z);
				
				PositionLineUnit = floor((FinalSample_theta-AngleOfFirstLine)/AngleStep) ;
				PositionRadialSampleUnit = floor((FinalSample_R-FirstSampleAxialPosition)/RadialStep) ;


				// test if the pixel is inside the sector

				if((PositionLineUnit>=0)&&(PositionRadialSampleUnit>=0)&&(PositionLineUnit<NTheta)&&(PositionRadialSampleUnit<NR))
				{
					LUTSIZE++;
				}
			}
		}
		
	// second pass
	
	else
	{

		// polar coordinates of output data (for pixel inside the sector)

		k = 0 ;
		for(int x=0;x<Nx;x++)
		{
			FinalSample_X = XAxis[x] ;

			for(int z=0;z<Nz;z++)
			{
				FinalSample_Z = ZAxis[z] ;
				FinalSample_R = sqrt(FinalSample_X*FinalSample_X + FinalSample_Z*FinalSample_Z);
				FinalSample_theta = atan2(FinalSample_X,FinalSample_Z);

				PositionLineUnit = floor((FinalSample_theta-AngleOfFirstLine)/AngleStep) ;
				PositionRadialSampleUnit = floor((FinalSample_R-FirstSampleAxialPosition)/RadialStep) ;


				// test if the pixel is inside the sector

				if((PositionLineUnit>=0)&&(PositionRadialSampleUnit>=0)&&(PositionLineUnit<NTheta)&&(PositionRadialSampleUnit<NR))
				{
					// fill lut
					LUT_coeffSum = 0 ;
					SeekOffSet = k*NbVoisTot ;
					LUT_X[k] = (short)(x) ;
					LUT_Z[k] = (short)(z) ;
					
					indexVoisin = 0 ;
					for(int t=0;t<2*NbVois;t++)
					{

						ClosiestLine =  PositionLineUnit + t - NbVois+1  ;
						ClosiestLine = MAX(ClosiestLine,0);
						ClosiestLine = MIN(ClosiestLine,NTheta-1);


						for(int r=0;r<2*NbVois;r++)
						{
							ClosiestAxialSample = PositionRadialSampleUnit + r - NbVois+1 ;
							ClosiestAxialSample = MAX(ClosiestAxialSample,0);
							ClosiestAxialSample = MIN(ClosiestAxialSample,NR-1);

							

							LUT_Vois1[indexVoisin + SeekOffSet] = (short)ClosiestLine ;
							LUT_Vois2[indexVoisin + SeekOffSet] = (short)ClosiestAxialSample ;
							Distance = (X0[ClosiestLine + ClosiestAxialSample*NTheta]-FinalSample_X) ;
							Distance *= Distance ;
							Distance +=  (Z0[ClosiestLine+ClosiestAxialSample*NTheta]-FinalSample_Z)*(Z0[ClosiestLine+ClosiestAxialSample*NTheta]-FinalSample_Z) ;
							
							LUT_coeff[indexVoisin + SeekOffSet] = exp(-Distance/Sigma) ;//1/(1+(Distance/Sigma));
							LUT_coeffSum += LUT_coeff[indexVoisin + SeekOffSet] ;
							indexVoisin++ ;

						}
					}

					for(int u=0 ;u<NbVoisTot;u++)
					{

						LUT_coeff[u+SeekOffSet] /=LUT_coeffSum ;
					}

					k++;
				}// endif

			}
		}
	}


	mexPrintf("LUTSIZE = %d\n",LUTSIZE);

	free(SINTHETA0);
	free(COSTHETA0);
	free(R0);
	free(Z0);
	free(X0);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


	int NR ;
	int NTheta ;
	float AngleOfFirstLine ;
	float AngleStep ;
	float FirstSampleAxialPosition  ;
	float RadialStep  ;
	float OffSet[2] ;
	float SizeOfImg[2];
	float OffSetX,OffSetZ, SizeOfImgX,SizeOfImgZ ;
	double * ptrLUTSIZE ;
	int Nx ;
	int Nz ; 
	int NbVois ;
	float SmoothDistance ;
	int LUTSIZE= 0 ;
	int InputLUTSIZE = 0 ;
	bool ComputeOffSetAndSize ;
	bool SecondPass = false ;
	int status ;
	int ScanconversionType ;

	double * parameters ;

	float * XAxis ;
	float * ZAxis ;

	float * LUT_coeff ;
	short * LUT_X ;
	short * LUT_Z ;
	short * LUT_Vois1 ;
	short * LUT_Vois2 ;

	int NbVoisTot ;

	parameters = (double*)mxGetData(prhs[0]) ;

	ScanconversionType = getIntValue(parameters,"ScanconversionType",scanconversionparameters,(int)parameters[1]); // not used yet
	SecondPass = getBoolValue(parameters,"SecondPass",scanconversionparameters,(int)parameters[1]);
	ComputeOffSetAndSize = getBoolValue(parameters,"ComputeOffSetAndSize",scanconversionparameters,(int)parameters[1]);
	NR = getIntValue(parameters,"NR",scanconversionparameters,(int)parameters[1]);
	NTheta = getIntValue(parameters,"NTheta",scanconversionparameters,(int)parameters[1]);
	AngleOfFirstLine = getFloatValue(parameters,"AngleOfFirstLine",scanconversionparameters,(int)parameters[1]);
	AngleStep = getFloatValue(parameters,"AngleStep",scanconversionparameters,(int)parameters[1]);
	FirstSampleAxialPosition = getFloatValue(parameters,"FirstSampleAxialPosition",scanconversionparameters,(int)parameters[1]);
	RadialStep = getFloatValue(parameters,"RadialStep",scanconversionparameters,(int)parameters[1]);
	OffSetX = getFloatValue(parameters,"OffSetX",scanconversionparameters,(int)parameters[1]);
	OffSetZ = getFloatValue(parameters,"OffSetZ",scanconversionparameters,(int)parameters[1]);
	SizeOfImgX = getFloatValue(parameters,"SizeOfImgX",scanconversionparameters,(int)parameters[1]);
	SizeOfImgZ = getFloatValue(parameters,"SizeOfImgZ",scanconversionparameters,(int)parameters[1]);
	Nx = getIntValue(parameters,"Nx",scanconversionparameters,(int)parameters[1]);
	Nz = getIntValue(parameters,"Nz",scanconversionparameters,(int)parameters[1]);
	NbVois = getIntValue(parameters,"NbVois",scanconversionparameters,(int)parameters[1]);
	SmoothDistance = getFloatValue(parameters,"SmoothDistance",scanconversionparameters,(int)parameters[1]);

	NbVoisTot = (2*NbVois) ; // Voisins
	NbVoisTot *= NbVoisTot ;

	if(SecondPass)
	{
		InputLUTSIZE = (*(int*)mxGetData(prhs[6]));

    // output
    plhs[0] = mxCreateDoubleScalar(InputLUTSIZE);

    int out_idx = 1;
    {
      mwSize ndim = 2;
      const mwSize dims[2] = { NbVoisTot * InputLUTSIZE, 1 };
      plhs[out_idx] = mxCreateNumericArray( ndim, dims, mxSINGLE_CLASS, mxREAL );
      LUT_coeff = (float*)mxGetData(plhs[out_idx]) ;
      out_idx++;
    }
    {
      mwSize ndim = 2;
      const mwSize dims[2] = { InputLUTSIZE, 1 };
      plhs[out_idx] = mxCreateNumericArray( ndim, dims, mxINT16_CLASS, mxREAL );
      LUT_X = (short*)mxGetData(plhs[out_idx]) ;
      out_idx++;
    }
    {
      mwSize ndim = 2;
      const mwSize dims[2] = { InputLUTSIZE, 1 };
      plhs[out_idx] = mxCreateNumericArray( ndim, dims, mxINT16_CLASS, mxREAL );
      LUT_Z = (short*)mxGetData(plhs[out_idx]) ;
      out_idx++;
    }
    {
      mwSize ndim = 2;
      const mwSize dims[2] = { NbVoisTot * InputLUTSIZE, 1 };
      plhs[out_idx] = mxCreateNumericArray( ndim, dims, mxINT16_CLASS, mxREAL );
      LUT_Vois1 = (short*)mxGetData(plhs[out_idx]) ;
      out_idx++;
    }
    {
      mwSize ndim = 2;
      const mwSize dims[2] = { NbVoisTot * InputLUTSIZE, 1 };
      plhs[out_idx] = mxCreateNumericArray( ndim, dims, mxINT16_CLASS, mxREAL );
      LUT_Vois2 = (short*)mxGetData(plhs[out_idx]) ;
      out_idx++;
    }
	}
	else
	{
    // Set output
    plhs[0] = mxCreateDoubleScalar(0);
    ptrLUTSIZE = mxGetPr(plhs[0]);

		LUTSIZE = 0 ;
	}
	
	XAxis = (float*)malloc(Nx*sizeof(float));
	ZAxis = (float*)malloc(Nz*sizeof(float));

	if(!ComputeOffSetAndSize)
	{
		OffSet[0] = OffSetX ;
		OffSet[1] = OffSetZ ;
		SizeOfImg[0] = SizeOfImgX ;
		SizeOfImg[1] = SizeOfImgZ ;
	}

	CompLut(NR,NTheta,AngleOfFirstLine,AngleStep,FirstSampleAxialPosition,RadialStep,OffSet,SizeOfImg,Nx,Nz,NbVois,SmoothDistance,XAxis,ZAxis, ComputeOffSetAndSize,LUT_coeff,LUT_X, LUT_Z ,LUT_Vois1, LUT_Vois2,LUTSIZE);

	if(ComputeOffSetAndSize)
	{
		SetValue(parameters,"OffSetX", scanconversionparameters, (int)parameters[1], (double)OffSet[0]);
		SetValue(parameters,"OffSetZ", scanconversionparameters, (int)parameters[1], (double)OffSet[1]);
		SetValue(parameters,"SizeOfImgX", scanconversionparameters, (int)parameters[1], (double)SizeOfImg[0]);
		SetValue(parameters,"SizeOfImgZ", scanconversionparameters, (int)parameters[1], (double)SizeOfImg[1]);
		OffSetX = OffSet[0];
		OffSetZ = OffSet[1];
		SizeOfImgX = SizeOfImg[0];
		SizeOfImgZ = SizeOfImg[1];
	}

	if(SecondPass)
	{
		ComputeOffSetAndSize = false ;

		if(InputLUTSIZE == LUTSIZE)
		{
			CompLut(NR,NTheta,AngleOfFirstLine,AngleStep,FirstSampleAxialPosition,RadialStep,OffSet,SizeOfImg,Nx,Nz,NbVois,SmoothDistance,XAxis,ZAxis, ComputeOffSetAndSize,
        LUT_coeff,LUT_X, LUT_Z ,LUT_Vois1, LUT_Vois2,LUTSIZE);
      mexPrintf("lutsize 2nd pass : %d\n", LUTSIZE);
		
			status = 0 ;
		}
		else
		{
			status = 1 ;
		}


		free(XAxis);
		free(ZAxis);
	}
	else // just output the lut size
	{
		ptrLUTSIZE[0] = (int)(LUTSIZE);
	}	
}
