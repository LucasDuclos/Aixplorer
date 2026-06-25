/*=================================================================
 *
 * 
 *
 * The calling syntax is:
 *
 *		FilterCoef = GenerateWallFilter()
 *
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 2010 SuperSonic Imagine
 *=================================================================*/
/* $Revision: 1.0 $ */
#include <math.h>
#include "mex.h"
#include        <fstream>
#include		<stdlib.h>
#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif
#define MAXENSEMBLELENGTH           64


struct S_WallFilter_IOP
{
    int    FilterType,
           FilterOrder,
           EnsembleLength;
    float  FilterCoef[MAXENSEMBLELENGTH*MAXENSEMBLELENGTH],
           FilterNumer[MAXENSEMBLELENGTH],
           FilterDenom[MAXENSEMBLELENGTH];
    float  FilterCutoff,
           StopBandAtten;
};




///////////////////////////////////////////////////////////////////////////////
// local functions
void        Identity_WallFilterCoefficientGeneration(S_WallFilter_IOP *WallFilter_IOP);
void        LSQ_WallFilterCoefficientGeneration(S_WallFilter_IOP *WallFilter_IOP);
void        Matrix_Multipl(double *A, double *B, double *C, int K, int L, int M);
void        Matrix_Inversion_GaussJordan(double *A, double *AI, int N);
void        IIR_WallFilterCoefficientGeneration(S_WallFilter_IOP *WallFilter_IOP);
void        Butterworth_LowPassAnalog(double *QuadrNumer, double *QuadrDenom, double X, int L);
void        Chebyshev2_LowPassAnalog(double *QuadrNumer, double *QuadrDenom, double X, int L);
void        LowPass_To_HighPass(double *QuadrNumer, double *QuadrDenom, double PreWarpedCutOff,
                                int Order, int HalfOrder);
void        Analog_To_Digital_BilinearTransform(double *QuadrNumer, double *QuadrDenom,
                                                double &CoefGain, int Order, int HalfOrder);
void        Quadratic_To_IIRFilterCoef(float *IirNumer, float *IirDenom, double *QuadrNumer, double *QuadrDenom,
                                       double CoefGain, int Order, int HalfOrder);
double      asinh(double x);
void        IIR_To_MatrixFilterCoef(float *IirNumer, float *IirDenom, float *FilterCoef,  int EL, int Order);
void        IIR_To_StateSpace(float *IirNumer, float *IirDenom, double *A,  double *B, double *C, int Order);
void        Matrix_Toeplitz(double *Col, double *Row, double *T, int N, int M);

///////////////////////////////////////////////////////////////////////////////
//
//	WallFilterCoefficientGeneration_v1
//  Stand-alone C program for generating wall-filter matrix coefficients
//	Thanasis Loupas, 17 September 2007
//
//
//	INPUTS
//
//	WallFilter_IOP->FilterType
//			Wall filter type. Currently, the following types are supported:
//					1: Butterworth
//					2: Chebyshev-II
//					3: Least-squares Regression
//                                      4: Identity
//	WallFilter_IOP->FilterOrder
//			IIR filter order (filter types 1 & 2) or least-squares polynomial degree (filter type 3).
//
//	WallFilter_IOP->FilterCutoff
//			IIR filter design cutoff (filter types 1 & 2).
//					Not relevant for filter type 3.
//
//	WallFilter_IOP->StopBandAtten
//			IIR filter stopband attenuation in dB.
//					Used only for filter type 2.
//					Hard-coded value = 10*log10(2) =3.01 db used for filter type 1
//					Irrelevant for filter type 3.
//
//	WallFilter_IOP->EnsembleLength
//			Input ensemble length
//
//
//	OUTPUTS
//
//	WallFilter_IOP->FilterCoef
//		2D array of EnsembleLength rows by EnsembleLength columns, containing the
//		projection-initialized matrix coefficients of the  wall filter designed according to the above specs
//			Floating-point precision.
//			Passband gain = 0 dB.
//			Skipping of filter matrix rows to minimize transient-response artifacts not implemented here
//
//
//	WARNING
//  No input argument checking.
//  Users are responsible for providing proper filter design parameters,
//  so that the achieved filter matches the design specifications.
//
//
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
void Identity_WallFilterCoefficientGeneration(S_WallFilter_IOP *WallFilter_IOP)
{
    int EL=WallFilter_IOP->EnsembleLength;
    WallFilter_IOP->FilterNumer[0] = 1.;
    WallFilter_IOP->FilterDenom[0] = 1.;
    for (int y=0; y<EL; y++)
    {
        for (int x=0; x<EL; x++)
        {
            WallFilter_IOP->FilterCoef[x*EL+y] = (y==x?1.:0);;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
void LSQ_WallFilterCoefficientGeneration(S_WallFilter_IOP *WallFilter_IOP)
{
	double *A=NULL, *ATr=NULL, *Buf1=NULL,  *Buf2=NULL, *Buf3=NULL, *Buf4=NULL;
	int Order1=0, EL=0, y=0, x=0, ayx=0;

	Order1=WallFilter_IOP->FilterOrder+1;
	EL=WallFilter_IOP->EnsembleLength;
	//TODO reallocate necessary
    A=new double[EL*Order1];
	//TODO reallocate necessary
    ATr=new double[Order1*EL];
	//TODO reallocate necessary
    Buf1=new double[Order1*Order1];
	//TODO reallocate necessary
    Buf2=new double[Order1*Order1];
	//TODO reallocate necessary
    Buf3=new double[EL*Order1];
	//TODO reallocate necessary
    Buf4=new double[EL*EL];
	for (y=0; y<EL; y++)
	{
		for (x=0; x<Order1; x++) {
			A[y*Order1+x]=pow(y+1., (double) x);
			ATr[x*EL+y]=A[y*Order1+x];
		}
	}
	Matrix_Multipl(ATr,A,Buf1,Order1,EL,Order1);
	Matrix_Inversion_GaussJordan(Buf1,Buf2,Order1);
	Matrix_Multipl(A,Buf2,Buf3,EL,Order1,Order1);
	Matrix_Multipl(Buf3,ATr,Buf4,EL,Order1,EL);
	for (y=0; y<EL; y++)
	{
		for (x=0; x<EL; x++)
		{
			ayx=y*EL+x;
			if (y==x)
			{
				WallFilter_IOP->FilterCoef[x*EL+y] =(float) (1.-Buf4[ayx]);
			}
			else
			{
				WallFilter_IOP->FilterCoef[x*EL+y] =(float) -Buf4[ayx];
			}
		}
	}
    
	delete[] A;
	delete[] ATr;
	delete[] Buf1;
	delete[] Buf2;
	delete[] Buf3;
	delete[] Buf4;
}

///////////////////////////////////////////////////////////////////////////////
void Matrix_Multipl(double *A, double *B, double *C, int K, int L, int M)
/*
	Matrix multiplication z=x*y
	x:	K by L
	y:	L by M
	z:	K by M
*/
{
	int k=0, m=0, l=0;
	double Sum=0.;

	for (k=0; k<K; k++)
	{
		for (m=0; m<M; m++)
		{
			Sum=0.;
			for (l=0; l<L; l++)
			{
				Sum += A[k*L+l]*B[l*M+m];
			}
			C[k*M+m]=Sum;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
void Matrix_Inversion_GaussJordan(double *A, double *AI, int N)
{
	int *IndxC=NULL, *IndxR=NULL, *Ipiv=NULL;
	int j=0, i=0, k=0, l=0, ll=0, aj=0, irow=0, icol=0;
	double MaxV=0., Tmp=0., PivInv=0.;

	for (j=0; j<N; j++)
	{
		aj=j*N;
		for (i=0; i<N; i++)
		{
			AI[aj+i]=A[aj+i];
		}
	}
	//TODO reallocate necessary
    IndxC=new int[N];
	//TODO reallocate necessary
    IndxR=new int[N];
	//TODO reallocate necessary
    Ipiv=new int[N];
	for (j=0; j<N; j++)
	{
		Ipiv[j]=0;
	}
	for (i=0; i<N; i++)
	{
		MaxV=0.;
		for (j=0; j<N; j++)
		{
			if (Ipiv[j] != 1)
			{
				for (k=0; k<N; k++)
				{
					if (Ipiv[k] == 0)
					{
						if (fabs(AI[j*N+k]) >= MaxV)
						{
							MaxV=fabs(AI[j*N+k]);
							irow=j;
							icol=k;
						}
					}
					else if (Ipiv[k] > 1)
					{
						printf("Matrix_Inversion_GaussJordan Error 1\n");
					}
				}
			}
		}
		++(Ipiv[icol]);
		if (irow != icol)
		{
			for (l=0; l<N; l++)
			{
				Tmp=AI[irow*N+l];
				AI[irow*N+l]=AI[icol*N+l];
				AI[icol*N+l]=Tmp;
			}
		}
		IndxR[i]=irow;
		IndxC[i]=icol;
		if (AI[icol*N+icol] == 0.0)
		{
			printf("Matrix_Inversion_GaussJordan Error2");
		}
		PivInv=1.0/AI[icol*N+icol];
		AI[icol*N+icol]=1.;
		for (l=0; l<N; l++)
		{
			AI[icol*N+l] *= PivInv;
		}
		for (ll=0; ll<N; ll++)
		{
			if (ll != icol)
			{
				Tmp=AI[ll*N+icol];
				AI[ll*N+icol]=0.;
				for (l=0; l<N; l++)
				{
					AI[ll*N+l] -= AI[icol*N+l]*Tmp;
				}
			}
		}
	}

	for (l=N-1; l>=0; l--)
	{
		if (IndxR[l] != IndxC[l])
		{
			for (k=0; k<N; k++)
			{
				Tmp=AI[k*N+IndxR[l]];
				AI[k*N+IndxR[l]]=AI[k*N+IndxC[l]];
				AI[k*N+IndxC[l]]=Tmp;
			}
		}
	}
	delete[] IndxR;
	delete[] IndxC;
	delete[] Ipiv;
}


///////////////////////////////////////////////////////////////////////////////
void IIR_WallFilterCoefficientGeneration(S_WallFilter_IOP *WallFilter_IOP)
{
	int Order=0, HalfOrder=0, Order1=0;
	double *QuadrNumer=NULL, *QuadrDenom=NULL;
	double PreWarpedCutOff=0., CoefGain=0.;

	Order=WallFilter_IOP->FilterOrder;
	HalfOrder=(int)(((float)Order + 1.0F) / 2.0F);
	//TODO reallocate necessary
    QuadrNumer=new double[3*HalfOrder];
	//TODO reallocate necessary
    QuadrDenom=new double[3*HalfOrder];
	if (WallFilter_IOP->FilterType==1)
	{
		// For Butterworth, set StopbandAtten to defualt value 3.01
		WallFilter_IOP->StopBandAtten=10.*log10(2.);
		Butterworth_LowPassAnalog(QuadrNumer,QuadrDenom,WallFilter_IOP->StopBandAtten,Order);
	}
	else if (WallFilter_IOP->FilterType==2)
	{
		Chebyshev2_LowPassAnalog(QuadrNumer,QuadrDenom,WallFilter_IOP->StopBandAtten,Order);
	}
	else
	{
		printf("FilterType %d is currently not supported",WallFilter_IOP->FilterType);
		printf("Using default filter type =1 (Butterworth) instead");
		// For Butterworth, set StopbandAtten to defualt value 3.01
		WallFilter_IOP->StopBandAtten=10.*log10(2.);
		Butterworth_LowPassAnalog(QuadrNumer,QuadrDenom,WallFilter_IOP->StopBandAtten,Order);
	}
	PreWarpedCutOff = 2*tan(M_PI*WallFilter_IOP->FilterCutoff);
	LowPass_To_HighPass(QuadrNumer,QuadrDenom,PreWarpedCutOff,Order,HalfOrder);
	CoefGain=1.;
	Analog_To_Digital_BilinearTransform(QuadrNumer,QuadrDenom,CoefGain,Order,HalfOrder);
	Order1=Order+1;
    
    Quadratic_To_IIRFilterCoef(WallFilter_IOP->FilterNumer,WallFilter_IOP->FilterDenom,QuadrNumer,QuadrDenom,CoefGain,Order,HalfOrder);
	IIR_To_MatrixFilterCoef(WallFilter_IOP->FilterNumer,WallFilter_IOP->FilterDenom,WallFilter_IOP->FilterCoef,WallFilter_IOP->EnsembleLength,Order);
    delete[] QuadrNumer;
	delete[] QuadrDenom;
}

///////////////////////////////////////////////////////////////////////////////
void Butterworth_LowPassAnalog(double *QuadrNumer, double *QuadrDenom, double X, int L)
{
	int     CntNum=0, CntDen=0, k=0;
	double  E=0., R=0., Theta=0., Sigma=0., Omega=0.;

	E = sqrt(pow(10.,X/10.)-1.);
	R = pow(E,-1./L);
	CntNum=0;
	CntDen=0;
	if(L % 2)
	{
		QuadrNumer[CntNum++]=0.0;
		QuadrNumer[CntNum++]=0.0;
		QuadrNumer[CntNum++]=R;
		QuadrDenom[CntDen++]=0.0;
		QuadrDenom[CntDen++]=1.0;
		QuadrDenom[CntDen++]=R;
	}
	for (k=0; k<L/2; k++)
	{
		Theta=M_PI*(2.*k+L+1.)/(2.*L);
		Sigma=R*cos(Theta);
		Omega=R*sin(Theta);
		QuadrNumer[CntNum++]=0.0;
		QuadrNumer[CntNum++]=0.0;
		QuadrNumer[CntNum++]=Sigma*Sigma+Omega*Omega;
		QuadrDenom[CntDen++]=1.0;
		QuadrDenom[CntDen++]=-2 * Sigma;
		QuadrDenom[CntDen++]=Sigma*Sigma+Omega*Omega;
	}
}

///////////////////////////////////////////////////////////////////////////////
void Chebyshev2_LowPassAnalog(double *QuadrNumer, double *QuadrDenom, double X, int L)
{
	int     CntNum=0, CntDen=0, k=0;
	double  E=0., D=0., Phi=0., Sigma=0., Omega=0., A=0.;

	E = sqrt(pow(10.,X/10.)-1.);
	D = asinh(E)/L;
	CntNum=0;
	CntDen=0;
	if(L % 2)
	{
		QuadrNumer[CntNum++]=0.0;
		QuadrNumer[CntNum++]=0.0;
		QuadrNumer[CntNum++]=1./sinh(D);
		QuadrDenom[CntDen++]=0.0;
		QuadrDenom[CntDen++]=1.0;
		QuadrDenom[CntDen++]=1./sinh(D);
	}
	for (k=0; k<L/2; k++)
	{
		Phi=M_PI*(2.*k+1.)/(2.*L);
		Sigma=-sinh(D)*sin(Phi);
		Omega=-cosh(D)*cos(Phi);
		A=Omega*Omega+Sigma*Sigma;
		QuadrNumer[CntNum++]=1.;
		QuadrNumer[CntNum++]=0.;
		QuadrNumer[CntNum++]=1./(cos(Phi)*cos(Phi));
		QuadrDenom[CntDen++]=1.;
		QuadrDenom[CntDen++]=-2.*Sigma/A;
		QuadrDenom[CntDen++]=(Omega*Omega)/(A*A)+(Sigma*Sigma)/(A*A);
	}
}


///////////////////////////////////////////////////////////////////////////////
void LowPass_To_HighPass(double *QuadrNumer, double *QuadrDenom, double PreWarpedCutOff, int Order, int HalfOrder)
{
	int m0=0, m=0, n=0;

	if(Order % 2)
	{
		QuadrNumer[2]=PreWarpedCutOff*QuadrNumer[1]/QuadrNumer[2];
		QuadrNumer[1]=1.;
		QuadrDenom[2]=PreWarpedCutOff/QuadrDenom[2];
		m0=1;
	}
	else
	{
		m0=0;
	}
	for(m=m0; m<HalfOrder; m++)
	{
		n=3*m;
		QuadrNumer[n+1]*=(PreWarpedCutOff/QuadrNumer[n+2]);
		QuadrNumer[n+2]=PreWarpedCutOff*PreWarpedCutOff*QuadrNumer[n]/QuadrNumer[n+2];
		QuadrNumer[n]=1.;
		QuadrDenom[n+1]*=(PreWarpedCutOff/QuadrDenom[n+2]);
		QuadrDenom[n+2]=PreWarpedCutOff*PreWarpedCutOff*QuadrDenom[n]/QuadrDenom[n+2];
		QuadrDenom[n]=1.;
	}
}


///////////////////////////////////////////////////////////////////////////////
void Analog_To_Digital_BilinearTransform(double *QuadrNumer, double *QuadrDenom,
										 double &CoefGain, int Order, int HalfOrder)
{
	int     i0=0, i=0, j=0;
	double  fs=0., fs2=0.,fs4=0. ,N0=0., N1=0., N2=0., D0=0. ,D1=0. ,D2=0.;

	fs=1.;
	fs2=2*fs;
	fs4=fs2*fs2;
	i0=0;
	if(Order % 2)
	{
		N0=QuadrNumer[2]+QuadrNumer[1]*fs2;
		N1=QuadrNumer[2]-QuadrNumer[1]*fs2;
		D0=QuadrDenom[2]+QuadrDenom[1]*fs2;
		D1=QuadrDenom[2]-QuadrDenom[1]*fs2;
		QuadrNumer[0]=1.;
		QuadrNumer[1]=N1/N0;
		QuadrNumer[2]=0.;
		QuadrDenom[0]=1.;
		QuadrDenom[1]=D1/D0;
		QuadrDenom[2]=0.;
		CoefGain*=(N0/D0);
		i0=1;
	}
	for(i=i0; i<HalfOrder; i++)
	{
		j=3*i;
		N0=QuadrNumer[j]*fs4+QuadrNumer[j+1]*fs2+QuadrNumer[j+2];
		N1=2*(QuadrNumer[j+2]-QuadrNumer[j]*fs4);
		N2=QuadrNumer[j]*fs4-QuadrNumer[j+1]*fs2+QuadrNumer[j+2];
		D0=QuadrDenom[j]*fs4+QuadrDenom[j+1]*fs2+QuadrDenom[j+2];
		D1=2*(QuadrDenom[j+2]-QuadrDenom[j]*fs4);
		D2=QuadrDenom[j]*fs4-QuadrDenom[j+1]*fs2+QuadrDenom[j+2];
		QuadrNumer[j]=1.;
		QuadrNumer[j+1]=N1/N0;
		QuadrNumer[j+2]=N2/N0;
		QuadrDenom[j]=1.;
		QuadrDenom[j+1]=D1/D0;
		QuadrDenom[j+2]=D2/D0;
		CoefGain*=(N0/D0);
	}
}


///////////////////////////////////////////////////////////////////////////////
void Quadratic_To_IIRFilterCoef(float *IirNumer, float *IirDenom, double *QuadrNumer, double *QuadrDenom,
						     double CoefGain, int Order, int HalfOrder)
{
	int Order1=0, N=0, M=0, m=0, i=0, j=0, k=0, n=0;
	float *TmpNumer=NULL, *TmpDenom=NULL, SumNumer=0, SumDenom=0;

	Order1=Order+1;
	//TODO reallocate necessary
    TmpNumer=new float[Order1];
	//TODO reallocate necessary
    TmpDenom=new float[Order1];

	N=2;
	if(Order % 2)
	{
		M=N-1;
	}
	else
	{
		M=N;
	}
	for (m=0; m<=M; m++)
	{
		TmpNumer[m]=QuadrNumer[m];
		TmpDenom[m]=QuadrDenom[m];
	}
	if (HalfOrder<2)
	{
		for (m=0; m<=M; m++)
		{
			IirNumer[m]=TmpNumer[m]*CoefGain;
			IirDenom[m]=TmpDenom[m];
		}
		delete[] TmpNumer;
		delete[] TmpDenom;
		return;
	}

	for (i=1; i<HalfOrder; i++)
	{
		j=3*i;
		for (k=0; k<=M+N; k++)
		{
			SumNumer=0.;
			SumDenom=0.;
			for (n=0; n<=N; n++)
			{
				m=k-n;
				if ( (m>=0) && (m<=M) )
				{
					SumNumer+=QuadrNumer[j+n]*TmpNumer[m];
					SumDenom+=QuadrDenom[j+n]*TmpDenom[m];
				}
			}
			IirNumer[k]=SumNumer;
			IirDenom[k]=SumDenom;
		}
		M+=N;
		for (m=0; m<=M; m++)
		{
			TmpNumer[m]=IirNumer[m];
			TmpDenom[m]=IirDenom[m];
		}
	}
	for (i=0; i<=Order; i++)
	{
		IirNumer[i]*=CoefGain;
	}
	delete[] TmpNumer;
	delete[] TmpDenom;
}

///////////////////////////////////////////////////////////////////////////////
double asinh(double x)
{
	return log(x+sqrt(x*x+1.));
}


///////////////////////////////////////////////////////////////////////////////
void IIR_To_MatrixFilterCoef(float *IirNumer, float *IirDenom, float *FilterCoef,  int EL, int Order)
// Matrix coefficients for projection-initialized filter specified
// by num[1,..,N] / den[1,...,N] transfer function.
// E= 		ensemble_length
// N=		filter_order+1
// B:		(N-1)x(N-1)
// C:		(N-1)x1
// A:		1x(N-1)
// F:		Ex(N-1)
// FS:		1x(N-1)
// row: 	1xE
// col:		1XE
// Bk:		(N-1)x(N-1)
// G:		ExE
// Ft:		(N-1)xE
// matcoef:	ExE
{
	int y=0, x=0, ay=0, r=0;
	double *A=NULL, *B=NULL, *BK=NULL, *C=NULL, *F=NULL, *FS=NULL, *FTr=NULL, *G=NULL;
	double *Row=NULL, *Col=NULL, *Buf1=NULL, *Buf2=NULL, Tmp[1]={0};

	//TODO reallocate necessary
    A=new double[Order];
	//TODO reallocate necessary
    B=new double[Order*Order];
	//TODO reallocate necessary
    C=new double[Order];
	//TODO reallocate necessary
    F=new double[EL*Order];
	//TODO reallocate necessary
    FS=new double[Order];

	IIR_To_StateSpace(IirNumer,IirDenom,B,C,A,Order);
	y=0;
	for (x=0; x<Order; x++)
	{
		ay=y*Order;
		F[ay+x]=A[x];
		FS[x]=A[x];
	}

	//TODO reallocate necessary
    Buf1=new double[Order];
	for (y=1; y<EL; y++)
	{
		ay=y*Order;
		Matrix_Multipl(FS,B,Buf1,1,Order,Order);
		for (x=0; x<Order; x++)
		{
			F[ay+x]=Buf1[x];
			FS[x]=Buf1[x];
		}
	}
	delete[] Buf1;

	//TODO reallocate necessary
    Row=new double[EL];
	//TODO reallocate necessary
    Col=new double[EL];
	//TODO reallocate necessary
    Row[0]=(double)IirNumer[0];
	for (y=1; y<EL; y++)
	{
		Row[y]=0;
	}
	//TODO reallocate necessary
    Col[0]=(double)IirNumer[0];
	//TODO reallocate necessary
    BK=new double[Order*Order];
	for (y=0; y<Order; y++)
	{
		ay=y*Order;
		for (x=0; x<Order; x++)
		{
			if (x==y)
			{
				BK[ay+x]=1.;
			}
			else
			{
				BK[ay+x]=0.;
			}
		}
	}

	//TODO reallocate necessary
    Buf1=new double[Order*Order];
	//TODO reallocate necessary
    Buf2=new double[Order];
	for (r=0; r<EL-1; r++)
	{
		if (r>0)
		{
			Matrix_Multipl(BK,B,Buf1,Order,Order,Order);
			for (y=0; y<Order; y++)
			{
				ay=y*Order;
				for (x=0; x<Order; x++)
				{
					BK[ay+x]=Buf1[ay+x];
				}
			}
		}
		Matrix_Multipl(A,BK,Buf2,1,Order,Order);
		Matrix_Multipl(Buf2,C,Tmp,1,Order,1);
		Col[r+1]=Tmp[0];
	}
	delete[] Buf1;
	delete[] Buf2;

	//TODO reallocate necessary
    G=new double[EL*EL];
	Matrix_Toeplitz(Col,Row,G,EL,EL);
	//TODO reallocate necessary
    FTr=new double[Order*EL];
	for (y=0; y<EL; y++)
	{
		for (x=0; x<Order; x++)
		{
			FTr[x*EL+y]=F[y*Order+x];
		}
	}
	//TODO reallocate necessary
    Buf1=new double[Order*Order];
	//TODO reallocate necessary
    Buf2=new double[Order*Order];
	Matrix_Multipl(FTr,F,Buf1,Order,EL,Order);
	Matrix_Inversion_GaussJordan(Buf1,Buf2,Order);
	delete[] Buf1;
	//TODO reallocate necessary
    Buf1=new double[EL*Order];
	Matrix_Multipl(F,Buf2,Buf1,EL,Order,Order);
	delete[] Buf2;
	//TODO reallocate necessary
    Buf2=new double[EL*EL];
	Matrix_Multipl(Buf1,FTr,Buf2,EL,Order,EL);
	delete[] Buf1;
	//TODO reallocate necessary
    Buf1=new double[EL*EL];
	Matrix_Multipl(Buf2,G,Buf1,EL,EL,EL);
	for (y=0; y<EL; y++)
	{
		ay=y*EL;
		for (x=0; x<EL; x++)
		{
			FilterCoef[x*EL+y]=(float) (G[ay+x]-Buf1[ay+x]);
		}
	}

	delete[] A;
	delete[] B;
	delete[] BK;
	delete[] C;
	delete[] F;
	delete[] FS;
	delete[] FTr;
	delete[] G;
	delete[] Row;
	delete[] Col;
	delete[] Buf1;
	delete[] Buf2;

}


///////////////////////////////////////////////////////////////////////////////
void IIR_To_StateSpace(float *IirNumer, float *IirDenom, double *A,  double *B, double *C, int Order)
/*
	Transfer function to state-space conversion.
	The NUM and DEN vectors must be of equal length N, and without leading zeros.
	N=	filter_order+1
	a:		(N-1)x(N-1)
	b:		(N-1)x1
	c:		1x(N-1)
*/
{
	int Order1=0, y=0, ay=0, x=0;
	float *Numer=NULL, *Denom=NULL;

	Order1=Order+1;

	//TODO reallocate necessary
    Numer=new float[Order1];
	//TODO reallocate necessary
    Denom=new float[Order1];
	for (x=0; x<Order1; x++)
	{
		Numer[x]=IirNumer[x]/IirDenom[0];
	}
	for (x=0; x<Order; x++)
	{
		Denom[x]=IirDenom[x+1]/IirDenom[0];
	}
	y=0;
	for (x=0; x<Order; x++)
	{
		ay=y*Order;
		A[ay+x]=(double)-Denom[x];
	}
	for (y=1; y<Order; y++)
	{
		ay=y*Order;
		for (x=0; x<Order; x++)
		{
			if (x==(y-1))
			{
				A[ay+x]=1.;
			}
			else
			{
				A[ay+x]=0.;
			}
		}
	}
	B[0]=1.;
	for (y=1; y<Order; y++)
	{
		B[y]=0.;
	}
	for (x=0; x<Order; x++)
	{
		C[x]=(double)(Numer[x+1]-Denom[x]*Numer[0]);
	}
	delete[] Numer;
	delete[] Denom;
}



///////////////////////////////////////////////////////////////////////////////
void Matrix_Toeplitz(double *Col, double *Row, double *T, int N, int M)
/*
	Non-symmetric Toeplitz matrix t with c[] as its first column
	and r[] as its first row.
	x:	1x(N+M-1)
	t:	NxM
*/
{
	int i=0, j=0;
	double *Buf=NULL;

	//TODO reallocate necessary
    Buf=new double[N+M-1];
	j=0;
	for (i=M-1; i>=1; i--)
	{
		Buf[j]=Row[i];
		j++;
	}
	for (i=0; i<N; i++)
	{
		Buf[j]=Col[i];
		j++;
	}
	for (i=0; i<N; i++)
	{
		for (j=0; j<M; j++)
		{
			T[i*M+j]=Buf[i-j+M-1];
		}
	}
	delete[] Buf;
}

int WallFilterCoefficientGeneration(S_WallFilter_IOP *WallFilter_IOP)
{
	int ret=0;

	if ( (WallFilter_IOP->FilterType==1) || (WallFilter_IOP->FilterType==2) )
	{
		IIR_WallFilterCoefficientGeneration(WallFilter_IOP);
		ret=1;
	}
	else if (WallFilter_IOP->FilterType==3)
	{
		LSQ_WallFilterCoefficientGeneration(WallFilter_IOP);
		ret=1;
	}
        else if (WallFilter_IOP->FilterType==4)
        {
            Identity_WallFilterCoefficientGeneration(WallFilter_IOP);
            ret=1;
        }
	else
	{
		ret=-1;
	}
	return ret;

}
static void computefunction(
       //OUTPUT
       float*          FilterCoef,

       //INPUT
       int     FilterType,
       int     FilterOrder,
       float   FilterCutoff,
       float   FilterStopBandAtten,
       int     EnsembleLength
       )

{
    for (int i=0;i<EnsembleLength;i++)
    {
        for (int j=0;j<EnsembleLength;j++)
        {
            int idx = i+j*EnsembleLength;
            FilterCoef[idx]=0.0;
            if (i==j)
            {
                FilterCoef[idx]= 1.0;
            }
        }
    }
    struct S_WallFilter_IOP CurrentWallFilterStruct;
    
   

    CurrentWallFilterStruct.FilterType = FilterType;
    CurrentWallFilterStruct.FilterOrder = FilterOrder;
    CurrentWallFilterStruct.EnsembleLength = EnsembleLength;
    CurrentWallFilterStruct.FilterCutoff = FilterCutoff;
    CurrentWallFilterStruct.StopBandAtten = FilterStopBandAtten;

    int Ret = WallFilterCoefficientGeneration(&CurrentWallFilterStruct);
    
    if (Ret)
    {
        for (int i=0;i<EnsembleLength;i++)
        {
            for (int j=0;j<EnsembleLength;j++)
            {
                int idx = i+j*EnsembleLength;
                FilterCoef[idx]=CurrentWallFilterStruct.FilterCoef[idx];
            }
        }
    }
    else
    {
        // do nothing
    }
    
    return;
}
 
void mexFunction( int nlhs, mxArray *plhs[],int nrhs, const mxArray*prhs[] )
     
{ 
       //OUTPUT
       float*          FilterCoef;

       //INPUT
       int     FilterType;
       int     FilterOrder;
       float           FilterCutoff;
       float           FilterStopBandAtten;
       int     EnsembleLength;

    
        /* Check for proper number of arguments */

        if (nlhs != 1) 
        { 
            mexErrMsgTxt("One output argument required."); 
        }  

        if (nrhs != 5) 
        { 
            mexErrMsgTxt("Five input arguments required."); 
        }  

        FilterType          = (int)mxGetScalar(prhs[0]);
        FilterOrder         = (int)mxGetScalar(prhs[1]);
        FilterCutoff        = (float)mxGetScalar(prhs[2]);
        FilterStopBandAtten = (float)mxGetScalar(prhs[3]);
        EnsembleLength      = (int)mxGetScalar(prhs[4]);
            
        mwSize DIM[2];
        DIM[0] = EnsembleLength;
        DIM[1] = EnsembleLength;
        plhs[0] = (mxArray*) mxCreateNumericArray(2,DIM, mxSINGLE_CLASS, mxREAL);
        FilterCoef = (float *)mxGetPr(plhs[0]);

        /* Do the actual computations in a subroutine */
        computefunction(FilterCoef,FilterType,FilterOrder,FilterCutoff,FilterStopBandAtten,EnsembleLength);
   
    return;
}

