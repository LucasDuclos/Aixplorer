///////////////////////////////////////////////////////
//                                                   //
//  Copyright © 2005-2008 by SuperSonic Imagine, SA  //
//         Confidential - All Right Reserved         //
//                                                   //
///////////////////////////////////////////////////////
/* =============================================================
* Reconstruct.cpp
V18
* introduce capability to reconstruct only partially for multithreading.
if ReconInfo[0]is an integer, then reconstruction starts from 0 and for NRecon= ReconInfo[0] reconstructions.
if ReconInfo[0]is not an integer, then :
    - the integer part of ReconInfo[0] is the cumulated number of reconstructions for all threads : MaxRec.
    - the fractional part of ReconInfo[0] is : (StartRec+NRecon/65536)/65536 where StartRec is the index of the first reconstruction for that thread, and NRecon is the number of reconstructions for that thread.
*modify the accumulation mode code (unrolling) to be compatible with the evolution of GCC (SSE shift right with immediate argument)
*the mex input is changed to allow multithreading
V17
* Same as V16
V16
*  Introduce reconstruction types 0210308, 0210504, 0220504
* Introduce accumulation for planewave operation
V14 2008 03 04
* Introduce output types 100,101,102,110
V13 2007 02 21
* Correct bug for indexed mode and recondelay larger than zorigin minus firstsample. Corrected modes :0010100, 0010308, 0010508, 0020504, 0020100 : Passage de Recolut de non signe en signe
V12 2008 02 13
* Implement mode 01210308, 01210504 and 01220504 for THI
V11 2007 09 25
*include reconstruction mode 100504
* change immediate mode. bocomes immediate time, but index coef
* suppress unused modes
V10 2007 09 05
*  Change piezo usage so that no more channels than necessary are used when a probe with less piezos than reconchannels is used.
V9 2007 04 19
*  Updated deltaprefetch values to take care of new Data format linked to new DMA
*  Reorganize data from 'Region' and 'ReconInfo' to 'AcqInfo' and 'ReconInfo'
V8 Correct test versus NSamples.
V7
*   Integration du parametre SyntAcq pour formatage des donnees en deux tirs par firing, avec donnees conjointes pour les deux tirs
*   Gestion d’un Offset memoire entre les donnees pour deux  canaux consecutifs
*   Implementation du mode grille rectangulaire. En utilisant le code OutputType = 10 les calculs et la sortie se font selon une grille rectangulaire. Si l’on est steere le PitchPix est en realite la projection verticale du pitchpix
V6
*   Passage de ReconDelay en µsecondes. Pour etre coherent avec TX.Delay
V5
*  Amelioration ponderation
    Modification de AntennaGain pour ne plus moyenner le gain d’antenne par paquet de 4 pixels. Une entree de AntennaGain par pixel... Moins d’artefact pour les faibles profondeurs.
    Ajout d’un parametre de gain dans AntennaGain pour faire varier le gain en fonction de la profondeur. Ceci est important en flat ou la contrbution d’un point est prise en compte inversement proportionnellement à la racine carrée de la profondeur. On multiplie donc en flat les pixels reconstruits proportionnellement à la racine carree de la profondeur avant de les diviser par le gain d’antenne. (proportionnel a la profondeur)
    à AntennaGain a donc trois entrées par pixel : la gain d’antenne pour les canaux a droite de l’ouverture, le gain d’antenne pour les canaux a gauche de l’ouverture, et le gain global en profondeur.
*  Troncature en mode OuputType 2
    Les shorts sont tronques en sortie à 65535
V4  Release le 06 10 2006.
    Correction de bugs V3 :
    o   Mauvais calcul de Lut pour faibles profondeurs, lorsque tous les canaux ne contribuent pas à une region donnée.
    o   Revue de l’algorithme de ponderation en considerant separement les contributions des canaux a gauche et a droite.
    o   Revue de la taille allouée à la lut pour eviter plantage à ouverture maximale
    o   Correction de l’apodisation
    Nouvelles fonctionnalites
    o   Implementation de TxAlongRx
    o   Gestion automatique debrayable du choix des modes
    o   Implementation du folding
    o   Implementation du mode 100%
    o   Implementation de lut Immediate pour petites Lut resultantes
V3
    C’est une evolution de V2, mais avec une RecoLut qui redevient absolue, et non pas incrementale (4 fois plus grande) :

    GlobalLut est reduite aux informations de ReconMux, ReconAbsc, et ReconDelay, un entier de 4 octets  pour chaque reconstruction
    Chaque entrée de RecoLut comporte un short pour coder l’apodisation, et un autre pour coder le temps.
    L’ouverture est symetrique par rapport à l’abscisse a l’origine de chaque ligne. Nouveau parametre RxAlongTx à 1
V2
    Reduction de la taille des LUT.
    GlobalLut est reduite aux informations de ReconMux, ReconAbsc, et ReconDelay
    RecoLut code une entree en delta de temps sur un octet
    L’ouverture est symetrique par rapport à l’abscisse a l’origine de chaque ligne
    Sortie en (I,Q) float, Intensite float, ou sqrt(Intensité) en short. Nouveau parametre OutputType = 0, 1, ou 2.
V1
    L’allocation mémoire pour le tableau de calculs intermediaires IQSum est reduite en taille d’un facteur 8 (bug).
    La structure PiezoChannel est remplacée par la structure ReconMux, moins generale, mais adaptée à la structure de Mux de Rubi V1 et Rubi V2. Pas de gestion de Mux externe avec cette version.
V0
    Algorithme de depart, livré le 15 09 2006

* =============================================================*/

#include "mex.h"
#include "reconstruct.h"
#include <math.h>
#include <string>
#include <emmintrin.h>
#include <xmmintrin.h>
#include <assert.h> //make care, it's the common assert to be compliant with all pope usages, not the RUBI assert!!!

/* The Reconstruct CPP function. */
#include <stdio.h>

#ifdef __GNUC__ 
  #ifdef __LP64__
      #define __ASMGCC64_OK
  #endif
#endif      
#ifndef __ASMGCC64_OK 
  #ifdef __GNUC__
      #warning This system is not compliant with the 64 bits pope assembler optimization.
      #warning    => a simple intrinsics version is building.
  #endif
#endif

#define NB_PIEZO_MAX 256
//        #define __ASMGCC64_OK
        
namespace sp
{
    void SteerOutput ( double * IQArray, int NPixLine, int NCols, int JZero, int DeltaJ, int * Offset );
    void SteerOutputFloat ( float * IQtemp, float * IQArray, int NPixLine, int NCols, int JZero, int DeltaJ, int * Offset , bool accum );
    //void * sp::reconstruct ( void * arg );
};
using namespace sp;
/* The Reconstruct CPP function. */
void * sp::reconstruct ( void * arg )
{
    EReconMethod ReconMethod = E_INTRINSICS ;
    //EReconMethod ReconMethod = E_ASM ;
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Set Different Lut Pointers
    int ProbeSize, AcqInfoSize, ReconInfoSize, GlobalLutSize, AntennaGain128Size, InterCoefSize, IQSumSize;
    int dataFormat;
//     unsigned short *RecoLut;
    short *RecoLut;
    unsigned int *GlobalLut;
    int  *RecoLutMoinsStart;
    double *Probe, *AcqInfo, *ReconInfo;
    __m128 *AntennaGain128;
    __m128i *InterCoef, *IQSum128i;

//     #define MaxSliceGain 128            //Gain multiplier for output

    short * Lut ;
    short *RcvDataShort ;
    float * IQPix ;
    float * IQPixImag ;
    float * IQbufferR = NULL ;
    float * IQbufferC = NULL;

    if(arg!=NULL)
    {
        Lut = ((ArgReconstruct*)arg)->Lut ;
        RcvDataShort = ((ArgReconstruct*)arg)->RcvDataShort ;
        IQPix = ((ArgReconstruct*)arg)->IQPix ;
        IQPixImag = ((ArgReconstruct*)arg)->IQPixImag ;
        IQbufferR = ((ArgReconstruct*)arg)->IQbuffer ;
        dataFormat = ((ArgReconstruct*)arg)->dataFormat ;
    }
    else
	return 0;
   

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /* Get Structures parameters */
    // Set first set of parameters
    ProbeSize = 3*4;
    AcqInfoSize = 11*4;
    ReconInfoSize = 26*4;
    Probe = ( double * ) &Lut[0];
    AcqInfo = ( double * ) &Lut[ProbeSize];
    ReconInfo = ( double * ) &Lut[ProbeSize + AcqInfoSize];

    int NPiezos, NLines, NPixLine, NSlices, NPixSlice, NRecon, StartRec, MaxRec, PlaneWave, Mode, BandWidth, OutputType;
    int FirstSample, NSamples, NChannels, MaxChannelsOneSide, PipeOffset, LutSizeBytes, TxAlongRx;
    long long Temp;
    int SyntAcq, ChannelOffset;
    double PiezoWidth, PiezoPitch, XOrigin, LinePitch, SteerAngle,DemodFrequ, SoundSpeed;
    double ZOrigin, PitchPix, BandInterp, FootInterp, PeakDelay, RXApodOne, RXApodZero;

    // 'Probe'
    NPiezos = ( int ) Probe[0]; // integer unit. number of piezos on probe
    PiezoWidth = Probe[1];      // mm unit. Piezo width.
    PiezoPitch = Probe[2];      // mm unit. Piezo pitch.

    // 'AcqInfo'
    SteerAngle = AcqInfo[0] ;    // radian unit. Steering angle.
    BandWidth = ( int ) AcqInfo[1] ;  // 200 100 50. (I,Q) acquisition bandwidth.
    // FrequDivider not needed here
    DemodFrequ = AcqInfo[3] ;    // MHz unit. Demodulation frequency.
    SyntAcq = ( int ) AcqInfo[4];    // Synthetic acquisition flag. 0 or 1.
    ZOrigin = AcqInfo[5] ;       // mm unit. depth of first reconstructed points.
    // Depth not needed here
    // InitialTime not needed here
    FirstSample = ( int ) AcqInfo[8] ;   // integer. Index of first sample sent to PC (useful for partial transmission).
    NSamples = ( int ) AcqInfo[9] ;   // integer. Number of samples per channel.
    // NAcqChannels not needed here

    // 'ReconInfo

    MaxRec = ( int ) ReconInfo[0];  // size of Global Lut.
    if ( MaxRec == ReconInfo[0] )   // case of full reconstruction
    {
        NRecon = MaxRec;
        StartRec = 0;
    }
    else                            // case of partial reconstruction for multithreading
    {
        Temp = ( long long ) ( ( ReconInfo[0]-double ( MaxRec ) ) *65536*65536 );
        StartRec = ( int ) ( Temp>>16 );
        NRecon = ( int ) ( Temp & 0x000000000000ffff );
    }
//     mexPrintf("NRecon = %d\n", NRecon);
//     mexPrintf("StartRec = %d\n", StartRec);
//     mexPrintf("MaxRec = %d\n", MaxRec);
    PlaneWave = ( int ) ReconInfo[1] ;  // 0 or 1. 1 if plane wave reconstruction, 0 if focalized reconstruction
    NLines = ( int ) ReconInfo[2] ;   // integer unit. Number of lines in AcqInfo. 1 for plane wave, >= 1 for focused wave.
    XOrigin = ReconInfo[3] ;       // PiezoPitch unit. Left line abscissa offset/left edge of reference piezo. usually 0 for plane Wave.
    LinePitch = ReconInfo[4] ;     // PiezoPitch unit. Horizontal distance between two lines. 1 for normal resolution, 0.5 for double resolution.
    PitchPix = ReconInfo[5] ;      // Lambda unit. Distance between two pixels on reconstructed line.
    NPixLine = ( int ) ReconInfo[6] ;  // integer unit. Number of reconstructed slices
    PeakDelay = ReconInfo[7] ;          // µsec unit . Delay between start of transmit pulse, and maximum of transmit power.
    PipeOffset = ( int ) ReconInfo[8] ; // integer. demod cycle unit . Time corresponding to first sample, assumed to be a 'I'.
    NChannels = ( int ) ReconInfo[9] ; // integer. Number of acquisition channels.
    MaxChannelsOneSide = ( int ) ReconInfo[10] ;   // Maximum number of channels contributing on one side of a region.
    ChannelOffset = ( int ) ReconInfo[11]; // Channel offset in receive buffer.
    TxAlongRx = ( int ) ReconInfo[12];  // 1 = Receive aperture is set around firing abscissa.
    BandInterp = ReconInfo[13] ;         // real [0 100]. Low pass filter for interpolation. 97 or 65
    FootInterp = ReconInfo[14] ;         // real . 100 * Foot of Hanning interpolation. 50 or 14
    RXApodOne = ReconInfo[15];          // max aperture value for full gain (Arditty coefficient)
    RXApodZero = ReconInfo[16];         // max aperture value for a contributing channel
    OutputType = ( int ) ReconInfo[17] ;   // Type of output : 0=(I,Q)   1= intensity
    NPixSlice = ( int ) ReconInfo[18];  // integer unit. Number of reconstructed pixels per slice.
    LutSizeBytes = ( int ) ReconInfo[19] ; // Lut size in bytes
    Mode = ( int ) ReconInfo[20] ;       // Will be recomputed later
    SoundSpeed = ReconInfo[21] ;    // m/sec unit. Speed of sound.
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Compute NSlices : number of slices
    NSlices = int ( ceil ( NPixLine/double ( NPixSlice ) ) );
    //Set sizes of structures
    if ( BandWidth == 200 ) InterCoefSize = 2*4*16*8;
    if ( BandWidth == 100 ) InterCoefSize = 2*16*8*8;
    IQSumSize=2*NSlices*NPixSlice*NLines*8;
    GlobalLutSize=MaxRec*2;
    AntennaGain128Size = 3* ( NSlices*NPixSlice+4 ) *2;  // 1 float per pixel of center line on each side of line, plus one float for depth gain
    // Set second set of pointers
    InterCoef = ( __m128i * ) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize];
    IQSum128i = ( __m128i * ) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize + InterCoefSize];
    int ASMLut[20];
    int *pASMLut = ASMLut;


//     GlobalLut = (unsigned int *) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize +InterCoefSize + IQSumSize];
//     AntennaGain128 = (float *) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize + InterCoefSize + IQSumSize + GlobalLutSize];
//     RecoLutMoinsStart = (int *) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize + InterCoefSize + IQSumSize + GlobalLutSize + AntennaGain128Size];

    AntennaGain128 = ( __m128 * ) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize +InterCoefSize + IQSumSize];
    GlobalLut = ( unsigned int * ) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize +InterCoefSize + IQSumSize + AntennaGain128Size];
    RecoLutMoinsStart = ( int * ) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize + InterCoefSize + IQSumSize + AntennaGain128Size + GlobalLutSize];
/*
    printf("InterCoefSize: %d\n", InterCoefSize);
    printf("IQSumSize: %d\n", IQSumSize);
    printf("AntennaGain128Size: %d\n", AntennaGain128Size);
    printf("GlobalLutSize: %d\n", GlobalLutSize);
*/
//     RecoLut = (unsigned short *) &Lut[((ProbeSize + AcqInfoSize + ReconInfoSize + InterCoefSize + IQSumSize + AntennaGain128Size + GlobalLutSize + 2)& 0xFFFFFFF8) + 0x00000008];
    RecoLut = ( short * ) &Lut[ ( ( ProbeSize + AcqInfoSize + ReconInfoSize + InterCoefSize + IQSumSize + AntennaGain128Size + GlobalLutSize + 2 ) & 0xFFFFFFF8 ) + 0x00000008];
    long recon_lut_size_nb = LutSizeBytes/2 - ( ( ( ProbeSize + AcqInfoSize + ReconInfoSize + InterCoefSize + IQSumSize + AntennaGain128Size + GlobalLutSize + 2 ) & 0xFFFFFFF8 ) + 0x00000008 );

    ////////////////////////////////////////////////////////////////////////////////////////////////

//     mexPrintf ("Reconstruction Mode : 0x%08X\n", Mode);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Reconstruction
    ////////////////////////////////////////////////////////////////////////////////////////////////

    short *ChannelStartPtr=0, *ChannelStartPtrRight=0, *ChannelStartPtrLeft=0, *IQPixshort=0;
    int i, j, k, l, k1, Delay32, Time32=0, Mux, Absc, NChannelsLeft, NChannelsRight, Channel, Slice, RecoLutInd, IndCoef;
    int ChannelNext, MinRightLeft, MaxRightLeft, ChannelLeft, ChannelRight,CentralChannel, AccumOffset=0;
    double XRegionCenter, Lambda;
    char *dataprefetch;
    __m64 *IQSum64;
    __m128i *data;
    __m128i xmma128i, xmmb128i, xmmc128i, xmmd128i, xmme128i;
    __m128 xmma128, xmmb128, xmmc128, xmmd128, xmme128, xmmf128,xmmg128, xmmh128, MaxLeftGain128, MaxRightGain128;
    int NRec, DeltaPrefetch, SliceStart, JZero, DeltaJ;
    int AccumPackets, PostDivideBy2, NumAccum, NSSEAcq, NSSEPacket;
    int AbscNext;
    int weightingNChannelsRight[NB_PIEZO_MAX], weightingNChannelsLeft[NB_PIEZO_MAX]; 
    
    // By default, the apertures weigths equal 1
    for (k=0; k<NB_PIEZO_MAX; k++)
    {
        weightingNChannelsRight[k] = 1;
        weightingNChannelsLeft[k]  = 1;
    }


    // Casting pointer to different types
    IQSum64     = ( __m64 * )     &IQSum128i[0];
    IQPixshort = ( short * ) &IQPix[0];
    
    Lambda=SoundSpeed/DemodFrequ/1000;                      // Compute Lambda in mm
    DeltaJ = int ( 65536*PitchPix*Lambda*tan ( SteerAngle ) / ( PiezoPitch*LinePitch ) +0.5 );     //DeltaJ for DeltaZ = 1
    
    // pointer to complex part of IQbuffer
    if(OutputType==13)
    {      
        if(DeltaJ == 0) // we won't go in SteerOutput, accumulate while store IQPix -> switch to OutputType 3
        {
            IQbufferR = IQPix ;
            IQbufferC = IQPixImag ;
            OutputType = 3 ;
        }
        else // we need to store IQ data before steering for accumulation
        {
            IQbufferC = &IQbufferR[NPixLine*NLines*NRecon];
        }
    }
    else // no accumulation or mode 3
    {
        IQbufferR = IQPix ;
        IQbufferC = IQPixImag ;
    }
        
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Reconstruction
    //
    // first perform data accumulation in planewave if needed
    if ( PlaneWave>1 )          // Case where we are in plane wave and there is a need for accumulation
    {
//        printf ( "accumulation\n" );
        //  First recover parameters packed inside planewave
        //mexPrintf("PlaneWaveR = %0x\n", PlaneWave);
        NumAccum        = ( PlaneWave & 0x0000FFFF );
        //mexPrintf("NumAccumR = %d\n", NumAccum);
        AccumPackets    = ( PlaneWave & 0x00FF0000 ) >> 16;
        //mexPrintf("AccumPacketsR %d\n", AccumPackets);
        PostDivideBy2   = ( PlaneWave & 0xFF000000 ) >> 24;
        //mexPrintf("PostDivideBy2R %d\n", PostDivideBy2);
        PlaneWave = 1;
        if ( PostDivideBy2 > 0 )
        {
            // then process accumulation
            for ( k=0; k<NChannels/ ( 1+SyntAcq ); k++ ) // loop on acquisition channels
            {
                NSSEAcq = NSamples* ( 1+SyntAcq ) /8;
                data = ( __m128i * ) &RcvDataShort[ ( ChannelOffset/2 ) *k];      // Start address of data for this acquisition channel
                // case of a single step accumulation

                if ( AccumPackets==1 )
                {
                    for ( i=0; i<NSSEAcq; i++ )                                 // loop on SSE samples in depth
                    {
                        for ( j=1; j<NumAccum; j++ )                            // loop on the number of samples to accumulate
                        {
                            data[i] = _mm_add_epi16 ( data[i], data[i+NSSEAcq*j] );  // accumulation
                        }   // end loop on accumulation
                        data[i]=_mm_srai_epi16 ( data[i], 1 );        // post reduce value and store back
                    } // end loop of samples in depth
                }
                // case were there are several accumulation packets
                else
                {
                    NSSEPacket = NSSEAcq*NumAccum/AccumPackets;
                    for ( i=0; i<NSSEAcq; i++ )                                 // loop on SSE samples in depth
                    {
                        xmmb128i= _mm_setzero_si128();          // initialize result for this SSE
                        for ( l=0; l<AccumPackets; l++ )        // loop on accumulation packets
                        {
                            xmma128i = data[i + NSSEPacket*l];
                            for ( j=1; j<NumAccum/AccumPackets; j++ )                            // loop on the number of samples to accumulate
                            {
                                xmma128i = _mm_add_epi16 ( xmma128i, data[i+NSSEPacket*l+NSSEAcq*j] );  // accumulation
                            }   // end loop on accumulation for one packet
                            xmmb128i= _mm_add_epi16 ( xmmb128i, _mm_srai_epi16 ( xmma128i, 1 ) );     // accumulate packet result divided by 2
                        } // end loop of the different packets
                        data[i]=_mm_srai_epi16 ( xmmb128i, 1 );
                    }    // end loop of SSEs in depth
                }
            }
        }
        else                    // case there is no postdivision by 2
        {
            // then process accumulation
            for ( k=0; k<NChannels/ ( 1+SyntAcq ); k++ ) // loop on acquisition channels
            {
                NSSEAcq = NSamples* ( 1+SyntAcq ) /8;
                data = ( __m128i * ) &RcvDataShort[ ( ChannelOffset/2 ) *k];      // Start address of data for this acquisition channel
                // case of a single step accumulation

                if ( AccumPackets==1 )
                {
                    for ( i=0; i<NSSEAcq; i++ )                                 // loop on SSE samples in depth
                    {
                        for ( j=1; j<NumAccum; j++ )                            // loop on the number of samples to accumulate
                        {
                            data[i] = _mm_add_epi16 ( data[i], data[i+NSSEAcq*j] );  // accumulation
                        }   // end loop on accumulation
                    } // end loop of samples in depth
                }
                // case were there are several accumulation packets
                else
                {
                    NSSEPacket = NSSEAcq*NumAccum/AccumPackets;
                    for ( i=0; i<NSSEAcq; i++ )                                 // loop on SSE samples in depth
                    {
                        xmmb128i= _mm_setzero_si128();          // initialize result for this SSE
                        for ( l=0; l<AccumPackets; l++ )        // loop on accumulation packets
                        {
                            xmma128i = data[i + NSSEPacket*l];
                            for ( j=1; j<NumAccum/AccumPackets; j++ )                            // loop on the number of samples to accumulate
                            {
                                xmma128i = _mm_add_epi16 ( xmma128i, data[i+NSSEPacket*l+NSSEAcq*j] );  // accumulation
                            }   // end loop on accumulation for one packet
                            xmmb128i= _mm_add_epi16 ( xmmb128i, _mm_srai_epi16 ( xmma128i, 1 ) );     // accumulate packet result divided by 2
                        } // end loop of the different packets
                    }    // end loop of SSEs in depth
                }
            }
        }

    }
    // Set CentralChannel = 1 when the region is centered around the reference Piezo
    XRegionCenter = XOrigin + ( NLines-1 ) *LinePitch/2;
    if ( XRegionCenter == 0.5 )
    {
        CentralChannel = 1;
    }
    else
    {
        CentralChannel = 0;
    }
    
    // For plane (wave + TxAlongRx=1) only
    if (PlaneWave==1 && TxAlongRx == 1)
    {
        // If the lines to reconstruct stand on the left of the first piezo
        // of the probe, just do nothing
        Absc = ( GlobalLut[StartRec]& 0x0000FF00 ) >> 8;  // Absc = 2nd byte
        if (Absc == 0)
        {
            NRec = StartRec;
            AbscNext  = ( GlobalLut[NRec+1]& 0x0000FF00 ) >> 8;
            while (AbscNext == 0 && NRec < StartRec + NRecon -1)
            {
                weightingNChannelsRight[NRec] = 0;
                Absc = AbscNext;
                NRec++ ;
                AbscNext  = ( GlobalLut[NRec+1]& 0x0000FF00 ) >> 8;
            }
            NRec-- ;
            weightingNChannelsRight[NRec] = 1;
        }
        // If the lines to reconstruct stand on the right of the last 
        // piezo of the probe, just do nothing
        Absc = ( GlobalLut[StartRec + NRecon -1]& 0x0000FF00 ) >> 8;  // Absc = 2nd byte
        if (Absc == (NPiezos-1))
        {
            NRec = StartRec + NRecon -1;
            AbscNext  = ( GlobalLut[NRec-1]& 0x0000FF00 ) >> 8;
            while (AbscNext == (NPiezos-1) && NRec>0)
            {
                weightingNChannelsLeft[NRec] = 0;
                weightingNChannelsRight[NRec] = 0;
                Absc = AbscNext;
                NRec-- ;
                AbscNext  = ( GlobalLut[NRec-1]& 0x0000FF00 ) >> 8;
            }
            NRec++;
            weightingNChannelsLeft[NRec] = 1;
        }
    }

    //loop on all reconstructions
    for ( NRec = StartRec; NRec < StartRec + NRecon; NRec++ )
    {
        // Initialize Structure for accumulating data
        for ( i=0; i<2*NLines*NSlices*NPixSlice; i++ )
        {
            IQSum128i[i] = _mm_setzero_si128();
        }
        Mux = GlobalLut[NRec]& 0x000000FF;          // Mux = 1st byte
        Absc = ( GlobalLut[NRec]& 0x0000FF00 ) >> 8;  // Absc = 2nd byte
        Delay32 = GlobalLut[NRec]>> 16;               // Delay = 3rd and 4th bytes
//         mexPrintf("Delay32 = %d\n", Delay32);
        NChannelsLeft=Absc-Mux;                     // Number of Receiving piezos left of region
        // Number of receiving piezos right of region
        if ( NPiezos < NChannels )
        {
            NChannelsRight =NPiezos-NChannelsLeft;      // case where probe has lower number of elements than number of recon channels
        }
        else
        {
            NChannelsRight =NChannels-NChannelsLeft;    // case where probe has number of elements egal or greater than number of recon channels
        }
        // Clipping to max channels used on each side, and compute max gain on each side
        if ( NChannelsLeft > ( MaxChannelsOneSide - CentralChannel ) ) NChannelsLeft = ( MaxChannelsOneSide- CentralChannel ); // Because the central channel is considered right
        NChannelsLeft  *= weightingNChannelsLeft[NRec];
        MaxLeftGain128 =  _mm_set1_ps ( float ( NChannelsLeft ) );
        if ( NChannelsRight > MaxChannelsOneSide ) NChannelsRight = MaxChannelsOneSide;
        NChannelsRight *= weightingNChannelsRight[NRec];
        MaxRightGain128 =  _mm_set1_ps ( float ( NChannelsRight ) );

        RecoLutInd = 0;
        // Compute first channel position
        Channel = Absc;
        if ( Channel >= NChannels ) Channel = Channel - NChannels;   // Take care of Mux

        pASMLut[0] = NSlices;                                   //0x10308 0x1110308
        pASMLut[1] = NLines;                                    //0x10308
        pASMLut[2] = ChannelOffset/2;                           //0x10308 0x1110308
        pASMLut[3] = NChannels/ ( 1+SyntAcq );                  //0x10308 0x1110308
        pASMLut[4] = NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave );   //0x10308 0x1110308
        pASMLut[5] = ( NChannels* ( 1+SyntAcq ) );              //0x10308 0x1110308
        pASMLut[6] = NSamples;                                  //0x10308 0x1110308
        pASMLut[7] = NChannelsRight;                            //0x10308 0x1110308
        pASMLut[8] = NChannelsRight+NChannelsLeft;              //0x10308 0x1110308
        pASMLut[9] = NChannels;                                 //0x10308 0x1110308
        pASMLut[10] = RecoLutMoinsStart[0];                     //0x10308
        pASMLut[12] = Absc;                                     //0x10308 0x1110308
        pASMLut[13] = NChannelsLeft+CentralChannel;             //        0x1110308
        pASMLut[14] = CentralChannel;                           //        0x1110308
        pASMLut[15] = NSamples*(1+SyntAcq);                     //                  0x1210504
        pASMLut[16] = 2*NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave );   //        0x1210504


        switch ( Mode )
        {
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1110308 :       // 3 pixels computed per slice with folding at 100% BandWidth and immediate Lut. Same code for 0x1120304
                #ifdef __ASMGCC64_OK
                if ( ReconMethod == E_ASM )
                {
                    // Determine number of folded channels
                    if ( ( NChannelsLeft+CentralChannel ) <= NChannelsRight )
                    {
                        MinRightLeft = NChannelsLeft+CentralChannel;           // Number of folded channels
                        MaxRightLeft = NChannelsRight;          // Max Channel number on one side
                    }
                    else
                    {
                        MaxRightLeft = NChannelsLeft+CentralChannel;
                        MinRightLeft = NChannelsRight;
                    }
                    ////////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////////
                    // Case where there is a central channel
                    if ( CentralChannel ==1 )
                    {
                        __asm__ __volatile__
                        (
                            ////////////////////////////////////////////////////////////////////////////
                            // Loop on contributing piezos
                            " xor   %%r15,              %%r15           \n" //r10=0
                            " mov   48(%[pASMLut]),     %%r15d          \n" // Channel = Absc

                            "9:                                             \n" //main loop on channels
                            " cmp   36(%[pASMLut]),     %%r15d          \n" //if Channel < NChannel,
                            " jb    10f                                 \n" //jump after the next instruction
                            " sub   36(%[pASMLut]),     %%r15d          \n" //else Channel = Channel - NChannels

                            "10:                                            \n"
                            //* Compute the following line
                            //* ChannelStartPtr= (short *) &RcvDataShort[(ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq))) +
                            //* NSamples*(1+SyntAcq)*NRec*(1-PlaneWave) +
                            //* NSamples*(Channel/(NChannels/(1+SyntAcq)))];

                            //(ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq)))
                            " xor   %%r13,              %%r13   \n" //r13=0
                            " movl  %%r15d,             %%eax   \n" //get Channel ==> low dividend part
                            " xor   %%edx,              %%edx   \n" //high dividend part is null
                            " divl 12(%[pASMLut])               \n" //channel % get pASMLut[3] => remainder in edx
                            " mov   %%edx,              %%eax   \n"
                            " mull  8(%[pASMLut])               \n" //multiply result by pASMLut[2] => result in eax

                            //NSamples*(1+SyntAcq)*NRec*(1-PlaneWave) ==> precomputed constant
                            " mov  16(%[pASMLut]),      %%r13d  \n" //get pASMLut[4]
                            " add   %%eax,              %%r13d  \n" //add pASMLut[4]+previous result ==> in r13

                            //NSamples*(Channel/(NChannels/(1+SyntAcq)))
                            " xor   %%edx,              %%edx   \n" //high dividend part is null
                            " movl  %%r15d,             %%eax   \n" //get Channel ==> low dividend part
                            " divl  12(%[pASMLut])              \n" //channel  get pASMLut[3] => quotient in eax
                            " mull  24(%[pASMLut])              \n" //multiply result by pASMLut[6] (NSamples) => result in eax
                            " add   %%eax,              %%r13d  \n" //both previous result ==> in r13

                            //apply the indirection to RcvDataShort to create ChannelStartPtr
                            " mov   %[RcvDataShort],    %[ChannelStartPtr]\n" //and create the final RcvDataShortPtr
                            " add   %%r13,              %[ChannelStartPtr]\n"
                            " add   %%r13,              %[ChannelStartPtr]\n" //twice, because it's shorts

                            ////////////////////////////////////////////////////////////////////////////
                            //loop on the Slices for this region
                            " xor   %%r13,              %%r13 \n" //r13=0 //the lines counter
                            " prefetch 64(%[RecoLut],   %%rbx, 2) \n" //prefetch the lut to be keep on L2 after use

                            //compute index values
                            //compute slice start for r14
                            " mov  (%[RecoLut], %%rbx, 2),     %%r14w \n" //get SliceStart (to replace after by SliceStart = RecoLut[RecoLutInd++];
                            " lea  (%%r14, %%r14),     %%r14  \n" //= SliceStart*2
                            " mov  (%[pASMLut]),       %%r8d  \n" //get NSlices
                            " lea  (%%r8, %%r8),       %%r8   \n" //= NSlices*2
                            " mov   %%r13,             %%r9   \n" //get j
                            " imul  %%r8,              %%r9   \n" //= j*2*NSlices
                            " add   %%r9,              %%r14  \n" //= SliceStart*2 + j*2*NSlices
                            " add   %%r8,              %%r9   \n" //= j*2*NSlices + NSlices*2
                            " mov   %%r14,             %%r8   \n" //get a slice_start copy
                            " sub   %%r8,              %%r9   \n" //= slice_start - slice_stop ==> nb_loop
                            " lea  (%%r14, %%r14,2),   %%r14  \n" //= (SliceStart*2 + j*2*NSlices)*3    ==> slice_start
                            " lea  (%%r14, %%r14),     %%r14  \n" //= (SliceStart*2 + j*2*NSlices)*3*2
                            " lea  (      , %%r14,8),  %%r14  \n" //= (SliceStart*2 + j*2*NSlices)*3*2*8 ==> offset for &IQSum128i[slice_start];
                            //" lea  (%%r9d, %%r9d,2),   %%r9d  \n" //= (j*2*NSlices + 2*NSlices)*3 ==> slice_stop

                            //validate values and start loop
                            " mov %%r9,                %%r15    \n"
                            " add %[RecoLutInd],       %%r15d   \n"
                            " cmp   $0,                %%r9     \n" //compare loop value with 0, to do not enter in loop if not necessary
                            " jle 2f                            \n" //immediat output case

                            " inc                      %[RecoLutInd]\n" //nb_loop_counter ++

                            "1:                                         \n"
                            //prefetch the lut
                            //load and compute corrects IndCoef and Time32
                            " movw  (%[RecoLut], %%rbx, 2), %%r9w  \n" //get Time 32 Lut Value
                            " inc    %[RecoLutInd]                 \n" //nb_loop_counter ++
                            " movw   (%[RecoLut], %%rbx, 2),  %%ax  \n" //get IndCoef Lut Value
                            " shl     $4,                    %%ax  \n" //InterCoef = InterCoef *4 Because InterCoef type is xmm128

//                             " prefetch    (%[InterCoef], %%rax) \n"

                            //get data[0]
                            " movdqu (%[ChannelStartPtr], %%r9, 2),    %%xmm0  \n" //load the source value from pSrc + Time32*sizeof(short)
                            //get data[1]
                            " movdqu 16(%[ChannelStartPtr], %%r9, 2),  %%xmm1  \n" //load the source value from pSrc + Time32*sizeof(short) + 16

                            " movdqa (%[InterCoef], %%rax),      %%xmm8  \n" //load InterCoef from    pRecolut + IndCoef
                            " movdqa 16(%[InterCoef], %%rax),    %%xmm9  \n" //load the InterCoef from pRecolut + IndCoef+16 (to get the next 128 block)

                            //Compute slice[0]:
                            " movdqa %%xmm8,  %%xmm7                    \n"
                            " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between source and lut
                            " paddd  (%[IQSum128i], %%r14),     %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm7,     (%[IQSum128i],%%r14)  \n" //copy result to dest

                            //Compute slice[1]:
                            " movdqa %%xmm9,  %%xmm7                    \n"
                            " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between source and lut
                            " paddd  16(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest[Slice+1]
                            " movdqa  %%xmm7,   16(%[IQSum128i],%%r14)  \n" //copy result to Dest[Slice+1]

                            //intermediate calculs
                            " shufps  $0x4E,    %%xmm1,         %%xmm0 \n" //shuffle 32 bits values with 4E= (1,0,3,2)

                            //Compute slice[2]:
                            " movdqa %%xmm8,  %%xmm7                    \n"
                            " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                            " paddd  32(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm7,   32(%[IQSum128i],%%r14)  \n" //copy result to dest

                            //Compute slice[3]:
                            " movdqa %%xmm9,  %%xmm7                    \n"
                            " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                            " paddd  48(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm7,   48(%[IQSum128i],%%r14)  \n" //copy result to dest

                            //Compute slice[4]:
                            " pmaddwd %%xmm1,                   %%xmm8  \n" //_mm_madd_epi16 between result of intermediate part and lut
                            " paddd  64(%[IQSum128i], %%r14),   %%xmm8  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm8,   64(%[IQSum128i],%%r14)  \n" //copy result to dest

                            //Compute slice[5]:
                            " pmaddwd %%xmm1,                   %%xmm9  \n" //_mm_madd_epi16 between result of intermediate part and lut
                            " paddd  80(%[IQSum128i], %%r14),   %%xmm9  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm9,   80(%[IQSum128i],%%r14)  \n" //copy result to dest

                            " inc                        %[RecoLutInd]  \n" //nb_loop_counter ++

                            " add     $96,                      %%r14   \n" //increase loop

                            " cmp     %%r15d,            %[RecoLutInd]  \n"
                            " jb       1b                               \n" //and branche it if non zero equal
                            "2:                                             \n" //on se casse...


                        :   [IQSum128i]         "+r" ( IQSum128i ),
                            [ChannelStartPtr]   "+r" ( ChannelStartPtr ),
                            [RecoLutInd]        "+b" ( RecoLutInd ),
                            [RecoLut]           "+r" ( RecoLut ),
                            [InterCoef]         "+r" ( InterCoef ),
                            [pASMLut]           "+r" ( pASMLut )
//                                [debug_val]         "+m" (debug_val),
//                                [debug_val2]        "+m" (debug_val2)
                                    :   [RcvDataShort]      "m" ( RcvDataShort ),
                                    [Delay32]           "m" ( Delay32 )

                                    : "%r15", "%r14", "%r13", "%r12", "%r10", "%r9", "%r8", "%rax", "%rdx",
                                    "%xmm0", "%xmm1", "%xmm2", "%xmm7", "%xmm8", "%xmm9"
                                );
                    }           // End loop on Channels

                    ////////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////////
                    // Loop on folded channels
                    if (  CentralChannel < MinRightLeft )
                    {
                        __asm__ __volatile__
                        (
                            ////////////////////////////////////////////////////////////////////////////
                            // Loop on contributing piezos
                            //Reset all used registers
                            " xor   %%r15,      %%r15  \n"
                            " xor   %%r14,      %%r14  \n"
                            " xor   %%r13,      %%r13  \n"
                            " xor   %%r10,      %%r10  \n"
                            " xor   %%r14,      %%r14  \n"
                            " xor   %%r9,       %%r9   \n"
                            " xor   %%r8,       %%r8   \n"
                            " xor   %%rdx,       %%rdx \n"
                            " xor   %%rax,       %%rax \n"
                            " mov   $0xFFFFFFFF, %%eax \n"
                            " and   %%rax, %%rbx       \n"
                            " xor   %%rax,       %%rax \n"

                            " mov   %[CentralChannel],    %%r10d        \n" // get the start counter (k)

                        "9:                                             \n" // main loop
                            " xor   %%r15,              %%r15           \n" // r15=0 ChannelLeft
                            " xor   %%r13,              %%r13           \n" // r13=0
                            " mov   %[Absc],            %%r15d          \n" // ChannelLeft = Absc
                            " mov   %%r15,              %%r8            \n" // ChannelRight = ChannelLeft

                            " add   %%r10,              %%r8            \n" // ChannelRight += k
                            " sub   %%r10,              %%r15           \n" // ChannelLeft -= main_counter (twice because it will be increase one after)
                            " add   56(%[pASMLut]),     %%r15d          \n" // ChannelLeft += CentralChannel
                            " dec                       %%r15           \n" // ChannelLeft --

                            " cmp   36(%[pASMLut]),     %%r15d          \n" //if ChannelLeft < NChannel,
                            " jb    10f                                 \n" //jump after the next instruction
                            " sub   36(%[pASMLut]),     %%r15d          \n" //else ChannelLeft = ChannelLeft - NChannels
                        "10:                                            \n"
                            " cmp   36(%[pASMLut]),     %%r8d           \n" //if ChannelRight < NChannel,
                            " jb    11f                                 \n" //jump after the next instruction
                            " sub   36(%[pASMLut]),     %%r8d           \n" //else ChannelRight = ChannelRight - NChannels
                        "11:                                            \n"
                            //* Compute the following line
                            //*     ChannelStartPtrLeft= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( ChannelLeft% ( NChannels/ ( 1+SyntAcq ) ) ) +
                            //*     NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                            //*     NSamples* ( ChannelLeft/ ( NChannels/ ( 1+SyntAcq ) ) ) ];

                            //(ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq)))
                            " xor   %%r13,              %%r13   \n" //r13=0
                            " movl  %%r15d,             %%eax   \n" //get Channel ==> low dividend part
                            " xor   %%edx,              %%edx   \n" //high dividend part is null
                            " divl 12(%[pASMLut])               \n" //channel % get pASMLut[3] => remainder in edx
                            " mov   %%edx,              %%eax   \n"
                            " mull  8(%[pASMLut])               \n" //multiply result by pASMLut[2] => result in eax

                            //NSamples*(1+SyntAcq)*NRec*(1-PlaneWave) ==> precomputed constant
                            " mov  16(%[pASMLut]),      %%r13d  \n" //get pASMLut[4]
                            " add   %%eax,              %%r13d  \n" //add pASMLut[4]+previous result ==> in r13

                            //NSamples*(Channel/(NChannels/(1+SyntAcq)))
                            " xor   %%edx,              %%edx   \n" //high dividend part is null
                            " movl  %%r15d,             %%eax   \n" //get Channel ==> low dividend part
                            " divl  12(%[pASMLut])              \n" //channel  get pASMLut[3] => quotient in eax
                            " mull  24(%[pASMLut])              \n" //multiply result by pASMLut[6] (NSamples) => result in eax
                            " add   %%eax,              %%r13d  \n" //both previous result ==> in r13

                            //apply the indirection to RcvDataShort to create ChannelStartPtr
                            " mov   %[RcvDataShort],    %[ChannelStartPtrLeft]\n" //and create the final RcvDataShortPtr
                            " add   %%r13,              %[ChannelStartPtrLeft]\n"
                            " add   %%r13,              %[ChannelStartPtrLeft]\n" //twice, because it's shorts

                            //* Compute the following line
                            //*     ChannelStartPtrRight= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( ChannelRight% ( NChannels/ ( 1+SyntAcq ) ) ) +
                            //*     NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                            //*     NSamples* ( ChannelRight/ ( NChannels/ ( 1+SyntAcq ) ) ) ];

                            //(ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq)))
                            " xor   %%r13,              %%r13   \n" //r13=0
                            " movl  %%r8d,             %%eax   \n" //get Channel ==> low dividend part
                            " xor   %%edx,              %%edx   \n" //high dividend part is null
                            " divl 12(%[pASMLut])               \n" //channel % get pASMLut[3] => remainder in edx
                            " mov   %%edx,              %%eax   \n"
                            " mull  8(%[pASMLut])               \n" //multiply result by pASMLut[2] => result in eax

                            //NSamples*(1+SyntAcq)*NRec*(1-PlaneWave) ==> precomputed constant
                            " mov  16(%[pASMLut]),      %%r13d  \n" //get pASMLut[4]
                            " add   %%eax,              %%r13d  \n" //add pASMLut[4]+previous result ==> in r13

                            //NSamples*(Channel/(NChannels/(1+SyntAcq)))
                            " xor   %%edx,              %%edx   \n" //high dividend part is null
                            " movl  %%r8d,             %%eax   \n" //get Channel ==> low dividend part
                            " divl  12(%[pASMLut])              \n" //channel  get pASMLut[3] => quotient in eax
                            " mull  24(%[pASMLut])              \n" //multiply result by pASMLut[6] (NSamples) => result in eax
                            " add   %%eax,              %%r13d  \n" //both previous result ==> in r13

                            //apply the indirection to RcvDataShort to create ChannelStartPtr
                            " mov   %[RcvDataShort],    %[ChannelStartPtrRight]\n" //and create the final RcvDataShortPtr
                            " add   %%r13,              %[ChannelStartPtrRight]\n"
                            " add   %%r13,              %[ChannelStartPtrRight]\n" //twice, because it's shorts

                            ////////////////////////////////////////////////////////////////////////////
                            //loop on the Slices for this region
                            " xor   %%r13,              %%r13 \n" //r13=0 //the lines counter
                            " prefetch 64(%[RecoLut],   %%rbx, 2) \n" //prefetch the lut to be keep on L2 after use

                            " xor   %%r14,      %%r14 \n" //r14=0
                            " xor   %%r8,      %%r8   \n" //r8=0

                            //compute index values
                            //compute slice start for r14
                            " mov  (%[RecoLut], %%rbx, 2),     %%r14w \n" //get SliceStart (to replace after by SliceStart = RecoLut[RecoLutInd++];
                            " lea  (%%r14, %%r14),     %%r14  \n" //= SliceStart*2
                            " mov  (%[pASMLut]),       %%r8d  \n" //get NSlices
                            " lea  (%%r8, %%r8),       %%r8   \n" //= NSlices*2
                            " mov   %%r13,             %%r9   \n" //get j
                            " imul  %%r8,              %%r9   \n" //= j*2*NSlices
                            " add   %%r9,              %%r14  \n" //= SliceStart*2 + j*2*NSlices
                            " add   %%r8,              %%r9   \n" //= j*2*NSlices + NSlices*2
                            " mov   %%r14,             %%r8   \n" //get a slice_start copy
                            " sub   %%r8,              %%r9   \n" //= slice_start - slice_stop ==> nb_loop
                            " lea  (%%r14, %%r14,2),   %%r14  \n" //= (SliceStart*2 + j*2*NSlices)*3    ==> slice_start
                            " lea  (%%r14, %%r14),     %%r14  \n" //= (SliceStart*2 + j*2*NSlices)*3*2
                            " lea  (      , %%r14,8),  %%r14  \n" //= (SliceStart*2 + j*2*NSlices)*3*2*8 ==> offset for &IQSum128i[slice_start];
                            //" lea  (%%r9d, %%r9d,2),   %%r9d  \n" //= (j*2*NSlices + 2*NSlices)*3 ==> slice_stop

                            //validate values and start loop
                            " mov %%r9,                %%r15    \n"
                            " add %[RecoLutInd],       %%r15d   \n"
                            " cmp   $0,                %%r9     \n" //compare loop value with 0, to do not enter in loop if not necessary
                            " jle 2f                            \n" //immediat output case

                            " inc                      %[RecoLutInd]\n" //nb_loop_counter ++

                        "1:                                         \n"
                            //prefetch the lut
                            //load and compute corrects IndCoef and Time32
                            " movw  (%[RecoLut], %%rbx, 2), %%r9w  \n" //get Time 32 Lut Value
                            " inc    %[RecoLutInd]                 \n" //nb_loop_counter ++
                            " movw   (%[RecoLut], %%rbx, 2),  %%ax  \n" //get IndCoef Lut Value
                            " shl     $4,                    %%ax  \n" //InterCoef = InterCoef *4 Because InterCoef type is xmm128

//                             " prefetch    (%[InterCoef], %%rax) \n"

                            //get data[0]
                            " movdqu (%[ChannelStartPtrLeft], %%r9, 2),    %%xmm0  \n" //load the source value from pSrc + Time32*sizeof(short)
                            " movdqu (%[ChannelStartPtrRight], %%r9, 2),    %%xmm2  \n" //load the source value from pSrc + Time32*sizeof(short)
                            //get data[1]
                            " movdqu 16(%[ChannelStartPtrLeft], %%r9, 2),  %%xmm1  \n" //load the source value from pSrc + Time32*sizeof(short) + 16
                            " movdqu 16(%[ChannelStartPtrRight], %%r9, 2),    %%xmm3  \n" //load the source value from pSrc + Time32*sizeof(short)

                            " paddw %%xmm2,                     %%xmm0   \n" //add folded channels
                            " paddw %%xmm3,                     %%xmm1   \n" //add folded channels n+1

                            " movdqa (%[InterCoef], %%rax),      %%xmm8  \n" //load InterCoef from    pRecolut + IndCoef
                            " movdqa 16(%[InterCoef], %%rax),    %%xmm9  \n" //load the InterCoef from pRecolut + IndCoef+16 (to get the next 128 block)

                            //Compute slice[0]:
                            " movdqa %%xmm8,  %%xmm7                    \n"
                            " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between source and lut
                            " paddd  (%[IQSum128i], %%r14),     %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm7,     (%[IQSum128i],%%r14)  \n" //copy result to dest

                            //Compute slice[1]:
                            " movdqa %%xmm9,  %%xmm7                    \n"
                            " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between source and lut
                            " paddd  16(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest[Slice+1]
                            " movdqa  %%xmm7,   16(%[IQSum128i],%%r14)  \n" //copy result to Dest[Slice+1]

                            //intermediate calculs
                            " shufps  $0x4E,    %%xmm1,         %%xmm0 \n" //shuffle 32 bits values with 4E= (1,0,3,2)

                            //Compute slice[2]:
                            " movdqa %%xmm8,  %%xmm7                    \n"
                            " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                            " paddd  32(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm7,   32(%[IQSum128i],%%r14)  \n" //copy result to dest

                            //Compute slice[3]:
                            " movdqa %%xmm9,  %%xmm7                    \n"
                            " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                            " paddd  48(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm7,   48(%[IQSum128i],%%r14)  \n" //copy result to dest

                            //Compute slice[4]:
                            " pmaddwd %%xmm1,                   %%xmm8  \n" //_mm_madd_epi16 between result of intermediate part and lut
                            " paddd  64(%[IQSum128i], %%r14),   %%xmm8  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm8,   64(%[IQSum128i],%%r14)  \n" //copy result to dest

                            //Compute slice[5]:
                            " pmaddwd %%xmm1,                   %%xmm9  \n" //_mm_madd_epi16 between result of intermediate part and lut
                            " paddd  80(%[IQSum128i], %%r14),   %%xmm9  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm9,   80(%[IQSum128i],%%r14)  \n" //copy result to dest

                            " inc                        %[RecoLutInd]  \n" //nb_loop_counter ++

                            " add     $96,                      %%r14   \n" //increase loop

                            " cmp     %%r15d,            %[RecoLutInd]  \n"
                            " jb       1b                               \n" //and branche it if non zero equal
                        "2:                                             \n" //on se casse...
                            " inc                       %%r10d          \n" //increment the main counter
                            " cmp     %[MinRightLeft],  %%r10d          \n" //and compare to the out value
                            " jne      9b                              \n" //if not the case jump to the start loop

                            :   [IQSum128i]         "+r" ( IQSum128i ),
                                [ChannelStartPtrLeft]  "+r" ( ChannelStartPtrLeft),
                                [ChannelStartPtrRight] "+r" ( ChannelStartPtrRight),
                                [RecoLutInd]        "+b" ( RecoLutInd ),
                                [RecoLut]           "+r" ( RecoLut ),
                                [InterCoef]         "+r" ( InterCoef ),
                                [pASMLut]           "+r" ( pASMLut )
//                                                                [debug_val]         "+m" (debug_val),
//                                                                [debug_val2]        "+m" (debug_val2)
                            :   [RcvDataShort]      "m" ( RcvDataShort ),
                                [Delay32]           "m" ( Delay32 ),
                                [MinRightLeft]      "m" ( MinRightLeft ),
                                [MaxRightLeft]      "m" ( MaxRightLeft ),
                                [CentralChannel]    "m" ( CentralChannel),
                                [Absc]              "m" ( Absc )

                            : "%r15", "%r14", "%r13", "%r10", "%r9", "%r8", "%rax", "%rdx",
                              "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm7", "%xmm8", "%xmm9"
                        );
                    }

                    ////////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////////
                    // Loop on unfolded channels
                    if ( MinRightLeft < MaxRightLeft )
                    {
                        __asm__ __volatile__
                        (
                            ////////////////////////////////////////////////////////////////////////////
                            // Loop on contributing piezos
                            " xor   %%r15,      %%r15 \n" //r15=0
                            " xor   %%r14,      %%r14 \n" //r14=0
                            " xor   %%r13,      %%r13 \n" //r14=0
                            " xor   %%r10,      %%r10           \n" //r10=0
                            " xor   %%r14,      %%r14 \n" //r14=0
                            " xor   %%r9,       %%r9   \n" //r9=0
                            " xor   %%r8,       %%r8   \n" //r8=0
                            " xor   %%rdx,       %%rdx   \n" //r8=0
                            " xor   %%rax,       %%rax   \n" //r8=0
                            " mov   $0xFFFFFFFF, %%eax \n "
                            " and   %%rax, %%rbx \n"
                            " xor   %%rax,       %%rax   \n" //r8=0

                            " mov   %[MinRightLeft],    %%r10d          \n" // get the start counter

                            "9:                                             \n" // main loop
                            " xor   %%r15,              %%r15           \n" //r15=0
                            " xor   %%r13,              %%r13           \n" //r15=0
                            " mov   %[Absc],            %%r15d          \n" // Channel = Absc
                            " mov   28(%[pASMLut]),     %%r13d          \n" // get NChannelsRight
                            " cmp   52(%[pASMLut]),     %%r13d          \n" // compare it with NChannelsLeft+ CentralChannel
                            " jg    12f                                 \n" // and jump to 12 if NChannelsRight > NChannelsLeft+ CentralChannel
                            " sub   %%r10,              %%r15           \n" //Channel -= main_counter
                            " sub   %%r10,              %%r15           \n" //Channel -= main_counter (twice because it will be increase one after)
                            " add   56(%[pASMLut]),     %%r15d          \n" //Channel += CentralChannel
                            " dec                       %%r15           \n" //Channel --
                            "12:                                            \n"
                            " add   %%r10,              %%r15           \n" //Channel += main_counter
                            " cmp   36(%[pASMLut]),     %%r15d          \n" //if Channel < NChannel,
                            " jb    10f                                 \n" //jump after the next instruction
                            " sub   36(%[pASMLut]),     %%r15d          \n" //else Channel = Channel - NChannels

                            "10:                                            \n"
                            //* Compute the following line
                            //* ChannelStartPtr= (short *) &RcvDataShort[(ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq))) +
                            //* NSamples*(1+SyntAcq)*NRec*(1-PlaneWave) +
                            //* NSamples*(Channel/(NChannels/(1+SyntAcq)))];

                            //(ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq)))
                            " xor   %%r13,              %%r13   \n" //r13=0
                            " movl  %%r15d,             %%eax   \n" //get Channel ==> low dividend part
                            " xor   %%edx,              %%edx   \n" //high dividend part is null
                            " divl 12(%[pASMLut])               \n" //channel % get pASMLut[3] => remainder in edx
                            " mov   %%edx,              %%eax   \n"
                            " mull  8(%[pASMLut])               \n" //multiply result by pASMLut[2] => result in eax

                            //NSamples*(1+SyntAcq)*NRec*(1-PlaneWave) ==> precomputed constant
                            " mov  16(%[pASMLut]),      %%r13d  \n" //get pASMLut[4]
                            " add   %%eax,              %%r13d  \n" //add pASMLut[4]+previous result ==> in r13

                            //NSamples*(Channel/(NChannels/(1+SyntAcq)))
                            " xor   %%edx,              %%edx   \n" //high dividend part is null
                            " movl  %%r15d,             %%eax   \n" //get Channel ==> low dividend part
                            " divl  12(%[pASMLut])              \n" //channel  get pASMLut[3] => quotient in eax
                            " mull  24(%[pASMLut])              \n" //multiply result by pASMLut[6] (NSamples) => result in eax
                            " add   %%eax,              %%r13d  \n" //both previous result ==> in r13

                            //apply the indirection to RcvDataShort to create ChannelStartPtr
                            " mov   %[RcvDataShort],    %[ChannelStartPtr]\n" //and create the final RcvDataShortPtr
                            " add   %%r13,              %[ChannelStartPtr]\n"
                            " add   %%r13,              %[ChannelStartPtr]\n" //twice, because it's shorts

                            ////////////////////////////////////////////////////////////////////////////
                            //loop on the Slices for this region
                            " xor   %%r13,              %%r13 \n" //r13=0 //the lines counter
                            " prefetch 64(%[RecoLut],   %%rbx, 2) \n" //prefetch the lut to be keep on L2 after use

                            " xor   %%r14,      %%r14 \n" //r14=0
                            " xor   %%r8,      %%r8   \n" //r8=0

                            //compute index values
                            //compute slice start for r14
                            " mov  (%[RecoLut], %%rbx, 2),     %%r14w \n" //get SliceStart (to replace after by SliceStart = RecoLut[RecoLutInd++];
                            " lea  (%%r14, %%r14),     %%r14  \n" //= SliceStart*2
                            " mov  (%[pASMLut]),       %%r8d  \n" //get NSlices
                            " lea  (%%r8, %%r8),       %%r8   \n" //= NSlices*2
                            " mov   %%r13,             %%r9   \n" //get j
                            " imul  %%r8,              %%r9   \n" //= j*2*NSlices
                            " add   %%r9,              %%r14  \n" //= SliceStart*2 + j*2*NSlices
                            " add   %%r8,              %%r9   \n" //= j*2*NSlices + NSlices*2
                            " mov   %%r14,             %%r8   \n" //get a slice_start copy
                            " sub   %%r8,              %%r9   \n" //= slice_start - slice_stop ==> nb_loop
                            " lea  (%%r14, %%r14,2),   %%r14  \n" //= (SliceStart*2 + j*2*NSlices)*3    ==> slice_start
                            " lea  (%%r14, %%r14),     %%r14  \n" //= (SliceStart*2 + j*2*NSlices)*3*2
                            " lea  (      , %%r14,8),  %%r14  \n" //= (SliceStart*2 + j*2*NSlices)*3*2*8 ==> offset for &IQSum128i[slice_start];
                            //" lea  (%%r9d, %%r9d,2),   %%r9d  \n" //= (j*2*NSlices + 2*NSlices)*3 ==> slice_stop

                            //validate values and start loop
                            " mov %%r9,                %%r15    \n"
                            " add %[RecoLutInd],       %%r15d   \n"
                            " cmp   $0,                %%r9     \n" //compare loop value with 0, to do not enter in loop if not necessary
                            " jle 2f                            \n" //immediat output case

                            " inc                      %[RecoLutInd]\n" //nb_loop_counter ++

                            "1:                                         \n"
                            //prefetch the lut
                            //load and compute corrects IndCoef and Time32
                            " movw  (%[RecoLut], %%rbx, 2), %%r9w  \n" //get Time 32 Lut Value
                            " inc    %[RecoLutInd]                 \n" //nb_loop_counter ++
                            " movw   (%[RecoLut], %%rbx, 2),  %%ax  \n" //get IndCoef Lut Value
                            " shl     $4,                    %%ax  \n" //InterCoef = InterCoef *4 Because InterCoef type is xmm128

//                             " prefetch    (%[InterCoef], %%rax) \n"

                            //get data[0]
                            " movdqu (%[ChannelStartPtr], %%r9, 2),    %%xmm0  \n" //load the source value from pSrc + Time32*sizeof(short)
                            //get data[1]
                            " movdqu 16(%[ChannelStartPtr], %%r9, 2),  %%xmm1  \n" //load the source value from pSrc + Time32*sizeof(short) + 16

                            " movdqa (%[InterCoef], %%rax),      %%xmm8  \n" //load InterCoef from    pRecolut + IndCoef
                            " movdqa 16(%[InterCoef], %%rax),    %%xmm9  \n" //load the InterCoef from pRecolut + IndCoef+16 (to get the next 128 block)

                            //Compute slice[0]:
                            " movdqa %%xmm8,  %%xmm7                    \n"
                            " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between source and lut
                            " paddd  (%[IQSum128i], %%r14),     %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm7,     (%[IQSum128i],%%r14)  \n" //copy result to dest

                            //Compute slice[1]:
                            " movdqa %%xmm9,  %%xmm7                    \n"
                            " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between source and lut
                            " paddd  16(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest[Slice+1]
                            " movdqa  %%xmm7,   16(%[IQSum128i],%%r14)  \n" //copy result to Dest[Slice+1]

                            //intermediate calculs
                            " shufps  $0x4E,    %%xmm1,         %%xmm0 \n" //shuffle 32 bits values with 4E= (1,0,3,2)

                            //Compute slice[2]:
                            " movdqa %%xmm8,  %%xmm7                    \n"
                            " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                            " paddd  32(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm7,   32(%[IQSum128i],%%r14)  \n" //copy result to dest

                            //Compute slice[3]:
                            " movdqa %%xmm9,  %%xmm7                    \n"
                            " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                            " paddd  48(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm7,   48(%[IQSum128i],%%r14)  \n" //copy result to dest

                            //Compute slice[4]:
                            " pmaddwd %%xmm1,                   %%xmm8  \n" //_mm_madd_epi16 between result of intermediate part and lut
                            " paddd  64(%[IQSum128i], %%r14),   %%xmm8  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm8,   64(%[IQSum128i],%%r14)  \n" //copy result to dest

                            //Compute slice[5]:
                            " pmaddwd %%xmm1,                   %%xmm9  \n" //_mm_madd_epi16 between result of intermediate part and lut
                            " paddd  80(%[IQSum128i], %%r14),   %%xmm9  \n" //_mm_add_epi32 between previous result and Dest
                            " movdqa  %%xmm9,   80(%[IQSum128i],%%r14)  \n" //copy result to dest

                            " inc                        %[RecoLutInd]  \n" //nb_loop_counter ++

                            " add     $96,                      %%r14   \n" //increase loop

                            " cmp     %%r15d,            %[RecoLutInd]  \n"
                            " jb       1b                               \n" //and branche it if non zero equal
                            "2:                                             \n" //on se casse...
                            " inc                       %%r10d          \n" //increment the main counter
                            " cmp     %[MaxRightLeft],  %%r10d          \n" //and compare to the out value
                            " jne      9b                              \n" //if not the case jump to the start loop

                        :   [IQSum128i]         "+r" ( IQSum128i ),
                            [ChannelStartPtr]   "+r" ( ChannelStartPtr ),
                            [RecoLutInd]        "+b" ( RecoLutInd ),
                            [RecoLut]           "+r" ( RecoLut ),
                            [InterCoef]         "+r" ( InterCoef ),
                            [pASMLut]           "+r" ( pASMLut ) /*,
                                [debug_val]         "+m" (debug_val),
                                [debug_val2]        "+m" (debug_val2)*/
                        :   [RcvDataShort]      "m" ( RcvDataShort ),
                            [Delay32]           "m" ( Delay32 ),
                            [MinRightLeft]      "m" ( MinRightLeft ),
                            [MaxRightLeft]      "m" ( MaxRightLeft ),
                            [Absc]              "m" ( Absc )

                        : "%r15", "%r14", "%r13", "%r10", "%r9", "%r8", "%rax", "%rdx",
                          "%xmm0", "%xmm1", "%xmm2", "%xmm7", "%xmm8", "%xmm9"
                        );
                    }

                } //end ASM
                else
                #endif //end #ifdef __ASMGCC64_OK
                {
                    // Determine number of folded channels
                    if ( ( NChannelsLeft+CentralChannel ) <= NChannelsRight )
                    {
                        MinRightLeft = NChannelsLeft+CentralChannel;           // Number of folded channels
                        MaxRightLeft = NChannelsRight;          // Max Channel number on one side
                    }
                    else
                    {
                        MaxRightLeft = NChannelsLeft+CentralChannel;
                        MinRightLeft = NChannelsRight;
                    }
                    ////////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////////
                    // Case where there is a central channel
                    if ( CentralChannel ==1 )
                    {
                        Channel= Absc;
                        // Acquisition Channel
                        if ( Channel >= NChannels ) Channel = Channel - NChannels;   // Acquisition Channel right
                        ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                         NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                         NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];

                        ////////////////////////////////////////////////////////////////////////////////
                        // Start the line
                        SliceStart =RecoLut[RecoLutInd++];
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart ) *3; Slice < ( 2*NSlices ) *3; Slice += 6 )
                        {
                            Time32 = RecoLut[RecoLutInd++];
                            IndCoef = RecoLut[RecoLutInd++];
                            xmma128i =_mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] );
                            IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            xmmc128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] );
                            xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                            xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                        }   // End loop on slices
                    }           // End loop on Channels
                    ////////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////////
                    // Loop on folded channels
                    for ( k=CentralChannel; k < MinRightLeft; k++ )
                    {
                        ChannelRight = Absc + k;
                        ChannelLeft = Absc - k -1+ CentralChannel;
                        // Acquisition Channel
                        if ( ChannelRight >= NChannels ) ChannelRight = ChannelRight - NChannels;   // Acquisition Channel right
                        if ( ChannelLeft >= NChannels ) ChannelLeft = ChannelLeft - NChannels;   // Acquisition Channel left
                        ChannelStartPtrRight= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( ChannelRight% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                              NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                              NSamples* ( ChannelRight/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                        ChannelStartPtrLeft= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( ChannelLeft% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                             NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                             NSamples* ( ChannelLeft/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                        /*                    DeltaPrefetch = (
                                            (ChannelOffset/2)*(ChannelNext%(NChannels/(1+SyntAcq))) +
                                            NSamples*(ChannelNext/(NChannels/(1+SyntAcq))) -
                                            (ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq))) -
                                            NSamples*(Channel/(NChannels/(1+SyntAcq)))
                                                            )/8; */
                        ////////////////////////////////////////////////////////////////////////////////
                        // Start line
                        SliceStart =RecoLut[RecoLutInd++];
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart ) *3; Slice < ( 2*NSlices ) *3; Slice += 6 )
                        {
                            Time32 = RecoLut[RecoLutInd++];
                            IndCoef = RecoLut[RecoLutInd++];
                            //                             dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32] + DeltaPrefetch);
                            //                             _mm_prefetch(dataprefetch,  _MM_HINT_NTA);
                            xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtrRight[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtrLeft[ ( Time32 ) + 0] ) );
                            IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            xmmc128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtrRight[ ( Time32 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtrLeft[ ( Time32 ) + 8] ) );
                            xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                            xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                        }   // End loop on slices
                    }           // End loop on Channels
                    ////////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////////
                    // Loop on unfolded channels
                    for ( k=MinRightLeft; k < MaxRightLeft; k++ )
                    {
                        if ( NChannelsRight > ( NChannelsLeft+ CentralChannel ) )
                        {
                            Channel = Absc + k;
                        }
                        else
                        {
                            Channel= Absc - k -1 + CentralChannel;
                        }
                        // Acquisition Channel
                        if ( Channel >= NChannels ) Channel = Channel - NChannels;   // Acquisition Channel right
                        ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                         NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                         NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                        ////////////////////////////////////////////////////////////////////////////////
                        // start line
                        SliceStart =RecoLut[RecoLutInd++];
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart ) *3; Slice < ( 2*NSlices ) *3; Slice += 6 )
                        {
                            Time32 = RecoLut[RecoLutInd++];
                            IndCoef = RecoLut[RecoLutInd++];
                            //                             dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32] + DeltaPrefetch);
                            //                             _mm_prefetch(dataprefetch,  _MM_HINT_NTA);
                            xmma128i =_mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] );
                            IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            xmmc128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] );
                            xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                            xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                        }   // End loop on slices
                    }           // End loop on Channels
                }
                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1010308 :       // 3 pixels computed per slice without folding at 100% BandWidth and immediate Lut. Same code for 0x1020304
                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }

                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart =RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ) *3; Slice < ( 2*NSlices + j*2*NSlices ) *3; Slice += 6 )
                        {
                            Time32 = RecoLut[RecoLutInd++];
                            IndCoef = RecoLut[RecoLutInd++];
//                             dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32] + DeltaPrefetch);
//                             _mm_prefetch(dataprefetch,  _MM_HINT_NTA);
                            xmma128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] );
                            IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            xmmc128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] );
                            xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                            xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1210308 :       // 3 pixels computed per slice without folding at 100% BandWidth with accumulate and immediate Lut. Same code for 0x1020304
                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }

                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     2*NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +     // multiply by 2 because there are two firings per recon
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    AccumOffset=NSamples* ( 1+SyntAcq );
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart =RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ) *3; Slice < ( 2*NSlices + j*2*NSlices ) *3; Slice += 6 )
                        {
                            Time32 = RecoLut[RecoLutInd++];
                            IndCoef = RecoLut[RecoLutInd++];
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32+AccumOffset] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + AccumOffset] ) );
                            IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            xmmc128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8+AccumOffset] ) );
                            xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                            xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1410308 :       // 3 pixels computed per slice without folding at 100% BandWidth with accumulate of 3 firings and immediate Lut. 
                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }

                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     3*NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +     // multiply by 3 because there are three firings per recon
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    AccumOffset=NSamples* ( 1+SyntAcq );
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart =RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ) *3; Slice < ( 2*NSlices + j*2*NSlices ) *3; Slice += 6 )
                        {
                            Time32 = RecoLut[RecoLutInd++];
                            IndCoef = RecoLut[RecoLutInd++];
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32+AccumOffset] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + AccumOffset] ) );
                            xmma128i =_mm_add_epi16 (xmma128i, _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 2*AccumOffset] ) );
                            IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            xmmc128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8+AccumOffset] ) );
                            xmmc128i =_mm_add_epi16 (xmmc128i, _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8 + 2*AccumOffset] ) );
                            xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                            xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1610308 :       // 3 pixels computed per slice without folding at 100% BandWidth with accumulate of 4 firings and immediate Lut. 
                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }

                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     4*NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +     // multiply by 4 because there are four firings per recon
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    AccumOffset=NSamples* ( 1+SyntAcq );
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart =RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ) *3; Slice < ( 2*NSlices + j*2*NSlices ) *3; Slice += 6 )
                        {
                            Time32 = RecoLut[RecoLutInd++];
                            IndCoef = RecoLut[RecoLutInd++];
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32+AccumOffset] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + AccumOffset] ) );
                            xmma128i =_mm_add_epi16 (xmma128i, _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 2*AccumOffset] ) );
                            xmma128i =_mm_add_epi16 (xmma128i, _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 3*AccumOffset] ) );
                            IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            xmmc128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8+AccumOffset] ) );
                            xmmc128i =_mm_add_epi16 (xmmc128i, _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8 + 2*AccumOffset] ) );
                            xmmc128i =_mm_add_epi16 (xmmc128i, _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8 + 3*AccumOffset] ) );
                            xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                            xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x10308 :       // 3 pixels computed per slice without folding at 100% BandWidth
				// TO DO : debug the assembly code for that mode                
				//if ( ReconMethod == E_ASM )
                //{
                    //////////////////////////////////////////////////////////////////////////////////
                    //// Loop on lines
                    //__asm__ __volatile__
                    //(
                        //////////////////////////////////////////////////////////////////////////////
                        //// Loop on contributing piezos
                        //" xor   %%r10,              %%r10           \n" //r10=0
                        //"9:                                             \n" //main loop on channels
                        //" xor   %%r14,              %%r14           \n" //r14=0
                        //" xor   %%r15,              %%r15           \n" //r15=0
                        //" mov   %%r10d,             %%r15d          \n" //get the channel counter
                        //" cmp   %%r15d,             28(%[pASMLut])  \n" //if (k == NChannelsRight)
                        //" cmove 40(%[pASMLut]),     %[RecoLutInd]   \n" //RecoLutInd = RecoLutMoinsStart[0]

                        //" inc   %%r15d                              \n" //k1 = k +1
                        //" cmp   %%r15d,             32(%[pASMLut])  \n" //if (k1 == NChannelsRight+NChannelsLeft)
                        //" jne   11f                                 \n"
                        //" dec   %%r15d                              \n" //k1 = k1 -1
                        //"11:                                            \n"
                        //" mov   48(%[pASMLut]),     %%r14d          \n"
                        //" cmp   28(%[pASMLut]),     %%r15d          \n" //if (k1 < NChannelsRight)
                        //" jb    12f                                 \n" //jump to 12
                        //" sub   %%r15d,             %%r14d          \n"
                        //" sub   %%r15d,             %%r14d          \n"
                        //" add   28(%[pASMLut]),     %%r14d          \n"
                        //" dec   %%r14d                              \n"
                        //"12:                                            \n"
                        //" add   %%r15d,             %%r14d          \n"

                        //" cmp   36(%[pASMLut]),     %%r14d          \n" //if ChannelNext < NChannel,
                        //" jb    10f                                 \n" //jump after the next instruction
                        //" sub   36(%[pASMLut]),     %%r14d          \n" // else ChannelNext = ChannelNext - NChannels
                        //"10:                                            \n"
                        //" xor   %%r13,          %%r13  \n" //r13=0
                        //" movl  %[Channel],     %%eax  \n" //get Channel ==> low dividend part
                        //" xor   %%edx,          %%edx  \n" //high dividend part is null
                        //" divl 12(%[pASMLut])          \n" //channel % get pASMLut[3] => remainder in edx
                        //" mov   %%edx,          %%eax  \n"
                        //" mull  8(%[pASMLut])          \n" //multiply result by pASMLut[2] => result in eax
                        //" mov  16(%[pASMLut]),  %%r13d \n" //get pASMLut[4]
                        //" add   %%eax,          %%r13d \n" //add pASMLut[4]+previous result ==> in r13

                        //" xor   %%edx,          %%edx  \n" //high dividend part is null
                        //" movl  %[Channel],     %%eax  \n" //get Channel ==> low dividend part
                        //" divl  12(%[pASMLut])         \n" //channel  get pASMLut[5] => quotient in eax
                        //" mull  24(%[pASMLut])         \n" //multiply result by pASMLut[6] (NSamples) => result in eax
                        //" add   %%eax,          %%r13d \n" //both previous result ==> in r13

                        //" mov   %[RcvDataShort], %[ChannelStartPtr]\n" //and create the final RcvDataShortPtr
                        //" add   %%r13, %[ChannelStartPtr]\n"
                        //" add   %%r13, %[ChannelStartPtr]\n" //twice, because it's shorts

                        //" mov   %%r14d,         %[Channel] \n"  //Channel = ChannelNext

                        //////////////////////////////////////////////////////////////////////////////
                        ////loop on the Slices for this region
                        //" xor   %%r13,      %%r13 \n" //r13=0 //the lines counter
                        //" prefetch 64(%[RecoLut], %%rbx, 2) \n" //prefetch the lut to be keep on L2 after use

                        //"3:\n " // line loop
                        //" xor   %%r15,      %%r15 \n" //r15=0
                        //" xor   %%r14,      %%r14 \n" //r14=0
                        //" xor   %%r9,      %%r9   \n" //r9=0
                        //" xor   %%r8,      %%r8   \n" //r8=0

                        ////compute index values
                        ////compute slice start for r14
                        //" mov  (%[RecoLut], %%rbx, 2),     %%r14w \n" //get SliceStart (to replace after by SliceStart = RecoLut[RecoLutInd++];
                        //" lea  (%%r14, %%r14),     %%r14 \n" //= SliceStart*2
                        //" mov  (%[pASMLut]),       %%r8d  \n" //get NSlices
                        //" lea  (%%r8, %%r8),       %%r8  \n" //= NSlices*2
                        //" mov   %%r13,             %%r9  \n" //get j
                        //" imul  %%r8,              %%r9   \n" //= j*2*NSlices
                        //" add   %%r9,              %%r14 \n" //= SliceStart*2 + j*2*NSlices
                        //" add   %%r8,              %%r9  \n" //= j*2*NSlices + NSlices*2
                        //" mov   %%r14,             %%r8  \n" //get a slice_start copy
                        //" sub   %%r8,              %%r9  \n" //= slice_start - slice_stop ==> nb_loop
                        //" lea  (%%r14, %%r14,2),   %%r14 \n" //= (SliceStart*2 + j*2*NSlices)*3    ==> slice_start
                        //" lea  (%%r14, %%r14),     %%r14 \n" //= (SliceStart*2 + j*2*NSlices)*3*2
                        //" lea  (      , %%r14,8),  %%r14 \n" //= (SliceStart*2 + j*2*NSlices)*3*2*8 ==> offset for &IQSum128i[slice_start];
                        ////" lea  (%%r9d, %%r9d,2),   %%r9d  \n" //= (j*2*NSlices + 2*NSlices)*3 ==> slice_stop

                        ////validate values and start loop
                        //" mov %%r9,                %%r15 \n"
                        //" add %[RecoLutInd],                  %%r15d  \n"
                        //" cmp   $0,                %%r9   \n" //compare loop value with 0, to do not enter in loop if not necessary
                        //" jle 2f                  \n" //immediat output case

                        //" inc                      %[RecoLutInd]   \n" //nb_loop_counter ++

                        //"1:\n"
                        ////prefetch the lut
                        ////load and compute corrects IndCoef and Time32
                        //" mov    %[Delay32],            %%r9d  \n" //Load Delay32
                        //" addw  (%[RecoLut], %%rbx, 2), %%r9w  \n" //Time 32 = Delay32 + Time32

                        //" mov       $0x1E,              %%rax  \n" //load 0x1E for IndCoef
                        //" and       %%r9,               %%rax  \n" //IndCoef = (Time32)& 0x1E

                        //" cmp       $7,                 %%rax  \n" //if IndCoef < $7
                        //" setle                         %%dl   \n" // prepare the value in rdx to decrease Time32

                        //" inc    %[RecoLutInd]                 \n" //nb_loop_counter ++
                        //" add   (%[RecoLut], %%rbx, 2),  %%ax  \n" //Load Lut IndCoef
                        //" shl     $4,                    %%ax  \n" //InterCoef = InterCoef *4 Because InterCoef type is xmm128

                        //" prefetch    (%[InterCoef], %%rax) \n"

                        ////shift Time 32 and mask it
                        //" shr    $4,               %%r9    \n" //Time32 = Time32 >> 4
                        //" and    $0xFFE,           %%r9    \n" //Time32 = Time32 & 0xFFE

                        ////and decrease it following the previous result in rdx
                        //" sub    %%rdx,             %%r9   \n" //if (Indcoef > 7) Time32 --

                        ////get data[0]
                        //" movdqu (%[ChannelStartPtr], %%r9, 2),    %%xmm0  \n" //load the source value from pSrc + Time32*sizeof(short)
                        ////get data[1]
                        //" movdqu 16(%[ChannelStartPtr], %%r9, 2),  %%xmm1  \n" //load the source value from pSrc + Time32*sizeof(short) + 16

                        //" movdqa (%[InterCoef], %%rax),      %%xmm8  \n" //load InterCoef from    pRecolut + IndCoef
                        //" movdqa 16(%[InterCoef], %%rax),    %%xmm9  \n" //load the InterCoef from pRecolut + IndCoef+16 (to get the next 128 block)

                        ////Compute slice[0]:
                        //" movdqa %%xmm8,  %%xmm7  \n"
                        //" pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between source and lut
                        //" paddd  (%[IQSum128i], %%r14),     %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        //" movdqa  %%xmm7,     (%[IQSum128i],%%r14)  \n" //copy result to dest

                        ////Compute slice[1]:
                        //" movdqa %%xmm9,  %%xmm7  \n"
                        //" pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between source and lut
                        //" paddd  16(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest[Slice+1]
                        //" movdqa  %%xmm7,   16(%[IQSum128i],%%r14)  \n" //copy result to Dest[Slice+1]

                        ////intermediate calculs
                        //" shufps  $0x4E,    %%xmm1,         %%xmm0 \n" //shuffle 32 bits values with 4E= (1,0,3,2)

                        ////Compute slice[2]:
                        //" movdqa %%xmm8,  %%xmm7  \n"
                        //" pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        //" paddd  32(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        //" movdqa  %%xmm7,   32(%[IQSum128i],%%r14)  \n" //copy result to dest

                        ////Compute slice[3]:
                        //" movdqa %%xmm9,  %%xmm7  \n"
                        //" pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        //" paddd  48(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        //" movdqa  %%xmm7,   48(%[IQSum128i],%%r14)  \n" //copy result to dest

                        ////Compute slice[4]:
                        //" pmaddwd %%xmm1,                   %%xmm8  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        //" paddd  64(%[IQSum128i], %%r14),   %%xmm8  \n" //_mm_add_epi32 between previous result and Dest
                        //" movdqa  %%xmm8,   64(%[IQSum128i],%%r14)  \n" //copy result to dest

                        ////Compute slice[5]:
                        //" pmaddwd %%xmm1,                   %%xmm9  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        //" paddd  80(%[IQSum128i], %%r14),   %%xmm9  \n" //_mm_add_epi32 between previous result and Dest
                        //" movdqa  %%xmm9,   80(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //" inc                        %[RecoLutInd]  \n" //nb_loop_counter ++

                        //" add     $96,                      %%r14   \n" //increase loop

                        //" cmp     %%r15d,            %[RecoLutInd]  \n"
                        //" jb       1b                               \n" //and branche it if non zero equal
                        //"2:                                             \n" //on se casse...
                        //" inc %%r13 \n " //increment line loop
                        //" cmp %%r13d, 4(%[pASMLut]) \n " //if < NLines
                        //" jb 3b \n " //jump to loop beginning*/

                        //" inc                       %%r10d          \n"
                        //" cmp       %%r10d,         32(%[pASMLut])  \n"
                        //" jne       9b                              \n"

                    //:   [IQSum128i]         "+r" ( IQSum128i ),
                        //[ChannelStartPtr]   "+r" ( ChannelStartPtr ),
                        //[RecoLutInd]        "+b" ( RecoLutInd ),
                        //[RecoLut]           "+r" ( RecoLut ),
                        //[InterCoef]         "+r" ( InterCoef ),
                        //[pASMLut]           "+r" ( pASMLut ),
                        //[Channel]           "+m" ( Channel ),
                        //[ChannelNext]       "+m" ( ChannelNext )
                                //:   [RcvDataShort]      "m" ( RcvDataShort ),
                                //[Delay32]           "m" ( Delay32 )

                                //: "%r15", "%r14", "%r13", "%r12", "%r10", "%r9", "%r8", "%rax", "%rdx",
                                //"%xmm0", "%xmm1", "%xmm2", "%xmm7", "%xmm8", "%xmm9"
                            //);
                //}
                //else
                {
                    for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )
                    {
                        // Loop on contributing piezos

                        if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                        // Compute next Channel for prefetch
                        k1 = k + 1;
                        if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                        if ( k1 < NChannelsRight )
                        {
                            ChannelNext = Absc + k1;
                        }
                        else
                        {
                            ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                        }

                        if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                        ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                         NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                         NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                        DeltaPrefetch = (
                                            ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                            NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                            ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                            NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                        ) /8;
                        Channel = ChannelNext;
                        ////////////////////////////////////////////////////////////////////////////////
                        // Loop on lines

                        for ( j=0; j<NLines; j++ )
                        {
                            SliceStart = RecoLut[RecoLutInd];

                            RecoLutInd++;
                            ////////////////////////////////////////////////////////////////////////////
                            //loop on the Slices for this region

                            for ( Slice = ( 2*SliceStart + j*2*NSlices ) *3; Slice < ( 2*NSlices + j*2*NSlices ) *3; Slice += 6 )
                            {
                                Time32 = RecoLut[RecoLutInd++];          // get time value in lambda*32
                                int Time32a = Time32 + Delay32;               // Time32 corrected by delay
                                IndCoef = ( Time32a & 0x001E );
                                if ( IndCoef >7 )           // Phase is >=4
                                {
                                    IndCoef += RecoLut[RecoLutInd++];
                                    Time32 = ( ( Time32a>>4 ) & 0xFFE );
                                    dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                                    _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                                    dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32+AccumOffset] + DeltaPrefetch );
                                    _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );

                                    xmma128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[Time32] ); //TEST
                                    IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                                    IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                                    xmmc128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] );
                                    xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                                    xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i ); //punpcklqdq
                                    IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                                    IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                                    IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                                    IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                                }
                                else                        // Phase is between 0 and 3
                                {

                                    Time32 = ( ( Time32a>>4 ) & 0xFFE ) -1;
                                    IndCoef += RecoLut[RecoLutInd++];
                                    xmma128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[Time32] ); //TEST
                                    IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                                    IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                                    xmmc128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] );
                                    xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                                    xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i ); //punpcklqdq
                                    IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                                    IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                                    IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                                    IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                                }
                            }   // End loop on slices
                        }
                    }       // End loop on Lines*/
                }           // End loop on Channels
                break;

                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x0210308 :       // 3 pixels computed per slice without folding at 100% BandWidth with accumulate
                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }

                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     2*NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +     // multiply by 2 because there is two firings per recon
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    AccumOffset=NSamples* ( 1+SyntAcq );
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart = RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ) *3; Slice < ( 2*NSlices + j*2*NSlices ) *3; Slice += 6 )
                        {
                            Time32 = RecoLut[RecoLutInd++];          // get time value in lambda*32
                            Time32 = Time32 + Delay32;               // Time32 corrected by delay
                            IndCoef = ( Time32 & 0x001E );//+ RecoLut[RecoLutInd++];
                            if ( IndCoef >7 )           // Phase is >=4
                            {
                                Time32 = ( Time32>>4 ) & 0xFFE;
                                IndCoef = IndCoef+ RecoLut[RecoLutInd++];
                                dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                                _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                                dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32+AccumOffset] + DeltaPrefetch );
                                _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                                //xmma128i = _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 0]);
                                xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + AccumOffset] ) );
                                IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                                IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                                //xmmc128i = _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 8]);
                                xmmc128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8+AccumOffset] ) );
                                xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                                xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                                IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                                IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                                IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                                IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                            }
                            else                        // Phase is between 0 and 3
                            {
                                Time32 = ( ( Time32>>4 ) & 0xFFE ) -1;
                                IndCoef = IndCoef+ RecoLut[RecoLutInd++];
//                                 dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32] + DeltaPrefetch);
//                                 _mm_prefetch(dataprefetch,  _MM_HINT_NTA);
//                                 dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32+AccumOffset] + DeltaPrefetch);
//                                 _mm_prefetch(dataprefetch,  _MM_HINT_NTA);
//                                 xmma128i = _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 0]);
                                xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + AccumOffset] ) );
                                IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                                IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
//                                 xmmc128i = _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 8]);
                                xmmc128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8+AccumOffset] ) );
                                xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                                xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                                IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                                IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                                IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                                IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                            }
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;

                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1010504 :       // 5 pixels computed per slice without folding at 100% BandWidth with immediate lut
                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }

                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart = RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ) *5; Slice < ( 2*NSlices + j*2*NSlices ) *5; Slice += 10 )
                        {
                            Time32 = RecoLut[RecoLutInd++];          // get time value in lambda*32
                            IndCoef = RecoLut[RecoLutInd++];
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            xmma128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] );
                            IQSum128i[Slice]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+0] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            xmme128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] );
                            xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                            xmmb128i=_mm_slli_si128 ( xmme128i, 12 );
                            xmmb128i=_mm_add_epi16 ( xmmb128i,xmma128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                            xmmc128i=_mm_slli_si128 ( xmme128i, 8 );
                            xmmc128i=_mm_add_epi16 ( xmmc128i,xmma128i );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                            xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                            xmmd128i=_mm_slli_si128 ( xmme128i, 4 );
                            xmmd128i=_mm_add_epi16 ( xmmd128i,xmma128i );
                            IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                            IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                            IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                            IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );

                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1210504 :       // 5 pixels computed per slice without folding at 100% BandWidth with accumulation and immediate lut
                #ifdef __ASMGCC64_OK
                if ( ReconMethod == E_ASM )
                {
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    __asm__ __volatile__
                    (
                        ////////////////////////////////////////////////////////////////////////////
                        // Loop on contributing piezos
                        " xor   %%r10,              %%r10           \n" //r10=0

                        "9:                                         \n" //main loop on channels
                        " xor   %%r14,              %%r14           \n" //r14=0
                        " xor   %%r15,              %%r15           \n" //r15=0
                        " mov   %%r10d,             %%r15d          \n" //get the channel counter
                        " cmp   %%r15d,             28(%[pASMLut])  \n" //if (k == NChannelsRight)
                        " cmove 40(%[pASMLut]),     %[RecoLutInd]   \n" //RecoLutInd = RecoLutMoinsStart[0]

                        " inc   %%r15d                              \n" //k1 = k +1
                        " cmp   %%r15d,             32(%[pASMLut])  \n" //if (k1 == NChannelsRight+NChannelsLeft)
                        " jne   11f                                 \n"
                        " dec   %%r15d                              \n" //k1 = k1 -1
                        "11:                                        \n"
                        " mov   48(%[pASMLut]),     %%r14d          \n"
                        " cmp   28(%[pASMLut]),     %%r15d          \n" //if (k1 < NChannelsRight)
                        " jb    12f                                 \n" //jump to 12
                        " sub   %%r15d,             %%r14d          \n"
                        " sub   %%r15d,             %%r14d          \n"
                        " add   28(%[pASMLut]),     %%r14d          \n"
                        " dec   %%r14d                              \n"
                        "12:                                        \n"
                        " add   %%r15d,             %%r14d          \n"

                        " cmp   36(%[pASMLut]),     %%r14d          \n" //if ChannelNext < NChannel,
                        " jb    10f                                 \n" //jump after the next instruction
                        " sub   36(%[pASMLut]),     %%r14d          \n" // else ChannelNext = ChannelNext - NChannels
                        "10:                                        \n"
                        " xor   %%r13,          %%r13  \n" //r13=0
                        " movl  %[Channel],     %%eax  \n" //get Channel ==> low dividend part
                        " xor   %%edx,          %%edx  \n" //high dividend part is null
                        " divl 12(%[pASMLut])          \n" //channel % get pASMLut[3] => remainder in edx
                        " mov   %%edx,          %%eax  \n"
                        " mull  8(%[pASMLut])          \n" //multiply result by pASMLut[2] => result in eax
                        " mov  64(%[pASMLut]),  %%r13d \n" //get pASMLut[16]
                        " add   %%eax,          %%r13d \n" //add pASMLut[16]+previous result ==> in r13

                        " xor   %%edx,          %%edx  \n" //high dividend part is null
                        " movl  %[Channel],     %%eax  \n" //get Channel ==> low dividend part
                        " divl  12(%[pASMLut])         \n" //channel  get pASMLut[5] => quotient in eax
                        " mull  24(%[pASMLut])         \n" //multiply result by pASMLut[6] (NSamples) => result in eax
                        " add   %%eax,          %%r13d \n" //both previous result ==> in r13

                        " mov   %[RcvDataShort], %[ChannelStartPtr]\n" //and create the final RcvDataShortPtr
                        " add   %%r13, %[ChannelStartPtr]\n"
                        " add   %%r13, %[ChannelStartPtr]\n" //twice, because it's shorts

                        " mov   %%r14d,         %[Channel] \n"  //Channel = ChannelNext
                        " xor   %%r13,      %%r13 \n" //r13=0 //the lines counter
                        " prefetch 64(%[RecoLut], %%rbx, 2) \n" //prefetch the lut to be keep on L2 after use

                        ////////////////////////////////////////////////////////////////////////////
                        "3:\n " // line loop
                        " xor   %%r14,      %%r14 \n" //r14=0
                        " xor   %%r8,      %%r8   \n" //r8=0

                        //compute index values
                        //compute slice start for r14
                        " mov  (%[RecoLut], %%rbx, 2),     %%r14w \n" //get SliceStart (to replace after by SliceStart = RecoLut[RecoLutInd++];

                        " lea  (%%r14, %%r14),     %%r14 \n" //= SliceStart*2
                        " mov  (%[pASMLut]),       %%r8d \n" //get NSlices
                        " lea  (%%r8, %%r8),       %%r8  \n" //= NSlices*2
                        " mov   %%r13,             %%r9  \n" //get j
                        " imul  %%r8,              %%r9  \n" //= j*2*NSlices
                        " add   %%r9,              %%r14 \n" //= SliceStart*2 + j*2*NSlices
                        " add   %%r8,              %%r9  \n" //= j*2*NSlices + NSlices*2
                        " mov   %%r14,             %%r8  \n" //get a slice_start copy
                        " sub   %%r8,              %%r9  \n" //= slice_start - slice_stop ==> nb_loop*2
                        " lea  (%%r14, %%r14,4),   %%r14 \n" //= (SliceStart*2 + j*2*NSlices)*5    ==> slice_start
                        " lea  (%%r14, %%r14),     %%r14 \n" //= (SliceStart*2 + j*2*NSlices)*5*2
                        " lea  (      , %%r14,8),  %%r14 \n" //= (SliceStart*2 + j*2*NSlices)*5*2*8 ==> offset for &IQSum128i[slice_start];
                        //" lea  (%%r9d, %%r9d,2),   %%r9d  \n" //= (j*2*NSlices + 2*NSlices)*3 ==> slice_stop

                        //validate values and start loop
                        " mov %%r9,                %%r15  \n"
                        " add %[RecoLutInd],       %%r15d \n"

                        " inc                 %[RecoLutInd]   \n" //nb_loop_counter ++

                        " cmp   $0,                %%r9   \n" //compare loop value with 0, to do not enter in loop if not necessary
                        " jle 2f                          \n" //immediat output case

                        ////////////////////////////////////////////////////////////////////////////
                        "1:\n"  //slice loop
                        " xor   %%r9,              %%r9           \n" //r9=0
                        " xor   %%rax,              %%rax           \n" //rax=0

                        //prefetch the lut
                        //load and compute corrects IndCoef and Time32
                        " movw  (%[RecoLut], %%rbx, 2),             %%r9w  \n" //get Time32

                        " inc    %[RecoLutInd]                          \n" //nb_loop_counter ++
                        " movw  (%[RecoLut], %%rbx, 2),             %%ax  \n" //get IndCoef
                        " shl     $4,                               %%ax  \n" //IndCoef = IndCoef <<4 Because InterCoef type is xmm128

                        //get InterCoef[0]
                        " movdqa (%[InterCoef], %%rax),             %%xmm8 \n" //load InterCoef from    pRecolut + IndCoef

                        //get data[0]
                        " movdqu (%[ChannelStartPtr], %%r9, 2),     %%xmm0 \n" //load the source value from pSrc + Time32*sizeof(short)
                        //get data[1]
                        " movdqu 16(%[ChannelStartPtr], %%r9, 2),   %%xmm1  \n" //load the source value from pSrc + Time32*sizeof(short) + 16

                         " add   60(%[pASMLut]),                     %%r9d  \n"  //compute accumOffset
                        //get data[0+accumOffset]
                         " movdqu (%[ChannelStartPtr], %%r9, 2),     %%xmm7 \n" //load the source value from pSrc + (Time32+AccumOffset)*sizeof(short)
                         " paddw  %%xmm7,                            %%xmm0 \n" //add both

                        //get data[1+accumOffset]
                         " movdqu 16(%[ChannelStartPtr], %%r9, 2),   %%xmm7 \n" //load the source value from pSrc + Time32*sizeof(short)
                         " paddw  %%xmm7,                            %%xmm1 \n" //add both

                        //Compute slice[0]:
                        " movdqa %%xmm8,                    %%xmm7  \n" //use positions dcba
                        " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between source and lut
                        " paddd  (%[IQSum128i], %%r14),     %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm7,     (%[IQSum128i],%%r14)  \n" //copy result to dest

                        //get InterCoef[1]
                       " movdqa 16(%[InterCoef], %%rax),   %%xmm9  \n" //load the InterCoef from pRecolut + IndCoef+16 (to get the next 128 block)

                        //Compute slice[1]:
                        " movdqa %%xmm9,                    %%xmm7  \n"
                        " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between source and lut
                        " paddd  16(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest[Slice+1]
                        " movdqa  %%xmm7,   16(%[IQSum128i],%%r14)  \n" //copy result to Dest[Slice+1]

                        //Compute slice[8]:
                        " movdqa  %%xmm8,                   %%xmm7  \n" // use positions hgfe
                        " pmaddwd %%xmm1,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  128(%[IQSum128i], %%r14),  %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm7,  128(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //Compute slice[9]:
                        " movdqa  %%xmm9,                   %%xmm7  \n"
                        " pmaddwd %%xmm1,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  144(%[IQSum128i], %%r14),  %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm7,  144(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //intermediate calculs
                        " pshufd $0x4E, %%xmm0,             %%xmm2  \n"
                        " punpcklqdq    %%xmm1,             %%xmm2  \n" //create the position for slides 4 & 5 => fedc

                        //Compute slice[4]:
                        " movdqa %%xmm8,                    %%xmm7  \n"
                        " pmaddwd %%xmm2,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  64(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm7,   64(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //Compute slice[5]:
                        " pmaddwd %%xmm9,                   %%xmm2  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  80(%[IQSum128i], %%r14),   %%xmm2  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm2,   80(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //intermediate calculs
                        " movss %%xmm1,                     %%xmm0 \n"  //create the position for slices 2 && 3 => edcb
                        " pshufd $0x39, %%xmm0,             %%xmm2 \n"

                        //Compute slice[2]:
                        " movdqa %%xmm8,                    %%xmm7  \n"
                        " pmaddwd %%xmm2,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  32(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm7,   32(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //Compute slice[3]:
                        " pmaddwd %%xmm9,                   %%xmm2  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  48(%[IQSum128i], %%r14),   %%xmm2  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm2,   48(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //intermediate calculs
                        " pshufd $0xFF, %%xmm0,             %%xmm0 \n"  //move the old MSB to LSB
                        " pslldq  $4,                       %%xmm1 \n"
                        " movss %%xmm0,                     %%xmm1 \n"  //create the position for slices 6 && 7 => gfed

                        //Compute slice[6]:
                        " pmaddwd %%xmm1,                   %%xmm8  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  96(%[IQSum128i], %%r14),   %%xmm8  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm8,   96(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //Compute slice[7]:
                        " pmaddwd %%xmm9,                   %%xmm1  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  112(%[IQSum128i], %%r14),  %%xmm1  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm1,  112(%[IQSum128i],%%r14)  \n" //copy result to dest*/

                        " inc                        %[RecoLutInd]  \n" //nb_loop_counter ++

                        " add     $160,                      %%r14   \n" //increase slice pointer (i128*10 = 160)

                        " cmp     %%r15d,            %[RecoLutInd]  \n" //nb_loop_counter = nb_loop
                        " jb       1b                               \n" //and branche it if non zero equal
                        "2:                                             \n" //on se casse...
                        " inc %%r13 \n " //increment line loop
                        " cmp  4(%[pASMLut]), %%r13d \n " //if < NLines
                        " jb 3b \n " //jump to loop beginning

                        " inc                       %%r10d          \n"
                        " cmp       %%r10d,         32(%[pASMLut])  \n"
                        " jne       9b                              \n"//*/
                        "30:\n"

                    :   [IQSum128i]         "+r" ( IQSum128i ),
                        [ChannelStartPtr]   "+r" ( ChannelStartPtr ),
                        [RecoLutInd]        "+b" ( RecoLutInd ),
                        [RecoLut]           "+r" ( RecoLut ),
                        [InterCoef]         "+r" ( InterCoef ),
                        [pASMLut]           "+r" ( pASMLut ),
                        [Channel]           "+m" ( Channel ),
                        [ChannelNext]       "+m" ( ChannelNext )
                        :   [RcvDataShort]      "m" ( RcvDataShort ),
                        [Delay32]           "m" ( Delay32 )

                    :   "%r15", "%r14", "%r13", "%r12", "%r10", "%r9", "%r8", "%rax", "%rdx",
                        "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm7", "%xmm8", "%xmm9"
                    );
//                     printf (" ASM debug 0x%08X \n", debug1);
                } //end ASM
                else
                #endif //end #ifdef __ASMGCC64_OK
                {
                    // Start loop on channels
                    for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                    {
                        if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                        // Compute next Channel for prefetch
                        k1 = k + 1;
                        if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                        if ( k1 < NChannelsRight )
                        {
                            ChannelNext = Absc + k1;
                        }
                        else
                        {
                            ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                        }

                        if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                        ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        2*NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +     // multiply by 2 because there are two firings per recon
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                        DeltaPrefetch = (
                                            ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                            NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                            ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                            NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                        ) /8;
                        Channel = ChannelNext;
                        AccumOffset=NSamples* ( 1+SyntAcq );
                        ////////////////////////////////////////////////////////////////////////////////
                        // Loop on lines
                        for ( j=0; j<NLines; j++ )
                        {
                            SliceStart = RecoLut[RecoLutInd++];

                            ////////////////////////////////////////////////////////////////////////////
                            //loop on the Slices for this region
                            for ( Slice = ( 2*SliceStart + j*2*NSlices ) *5; Slice < ( 2*NSlices + j*2*NSlices ) *5; Slice += 10 )
                            {
                                Time32 = RecoLut[RecoLutInd++];          // get time value in lambda*32
                                IndCoef = RecoLut[RecoLutInd++];
                                dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                                _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                                dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32+AccumOffset] + DeltaPrefetch );
                                _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                               xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 )] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 )+ AccumOffset] ) );
//                                 xmma128i =_mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] );
                                IQSum128i[Slice]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+0] );
                                IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                                xmme128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) +8+ AccumOffset] ) );
                                xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                xmmb128i=_mm_slli_si128 ( xmme128i, 12 );
                                xmmb128i=_mm_add_epi16 ( xmmb128i,xmma128i );
                                IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                                IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                                xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                xmmc128i=_mm_slli_si128 ( xmme128i, 8 );
                                xmmc128i=_mm_add_epi16 ( xmmc128i,xmma128i );
                                IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                                IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                                xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                xmmd128i=_mm_slli_si128 ( xmme128i, 4 );
                                xmmd128i=_mm_add_epi16 ( xmmd128i,xmma128i );
                                IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                                IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                                IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                                IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] ); //*/

                            }   // End loop on slices

                        }       // End loop on Lines
                    }           // End loop on Channels
                    break;


                }//end intrinsics
                break;

                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1410504 :       // 5 pixels computed per slice without folding at 100% BandWidth with accumulation of 3 firings and immediate lut
                // Start loop on channels
                for (k=0; k < NChannelsRight+NChannelsLeft; k++) {      // Loop on contributing piezos
                    if (k == NChannelsRight) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if (k1 == NChannelsRight+NChannelsLeft) k1 = k1 -1; // last Channel case
                    if (k1 < NChannelsRight) {
                        ChannelNext = Absc + k1;
                    }
                    else {
                        ChannelNext = Absc - (k1-NChannelsRight) -1;
                    }
                    
                    if (ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;    // Acquisition Channel
                    ChannelStartPtr= (short *) &RcvDataShort[(ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq))) +
                                                               3*NSamples*(1+SyntAcq)*NRec*(1-PlaneWave) +     // multiply by 3 because there are Three firings per recon
                                                                NSamples*(Channel/(NChannels/(1+SyntAcq)))];
                    DeltaPrefetch = (
                                    (ChannelOffset/2)*(ChannelNext%(NChannels/(1+SyntAcq))) +       
                                    NSamples*(ChannelNext/(NChannels/(1+SyntAcq))) -
                                    (ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq))) -
                                    NSamples*(Channel/(NChannels/(1+SyntAcq)))
                                    )/8;
                    Channel = ChannelNext;
                    AccumOffset=NSamples*(1+SyntAcq);
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for (j=0; j<NLines; j++) {
                        SliceStart = RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for (Slice = (2*SliceStart + j*2*NSlices)*5; Slice <(2*NSlices + j*2*NSlices)*5; Slice += 10) {
                            Time32 = RecoLut[RecoLutInd++];          // get time value in lambda*32
                            IndCoef = RecoLut[RecoLutInd++];
                            dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32] + DeltaPrefetch);
                            _mm_prefetch(dataprefetch,  _MM_HINT_NTA);
                            dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32+AccumOffset] + DeltaPrefetch);
                            _mm_prefetch(dataprefetch,  _MM_HINT_NTA);
                            xmma128i =_mm_add_epi16(_mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 0]), _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ AccumOffset]));
                            xmma128i =_mm_add_epi16(xmma128i, _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 2*AccumOffset]));
                            IQSum128i[Slice]=_mm_add_epi32( _mm_madd_epi16(xmma128i,InterCoef[IndCoef  ]), IQSum128i[Slice+0]);
                            IQSum128i[Slice+1]=_mm_add_epi32( _mm_madd_epi16(xmma128i,InterCoef[IndCoef+1]), IQSum128i[Slice+1]);
                            xmme128i =_mm_add_epi16(_mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 8]), _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+8+ AccumOffset]));
                            xmme128i =_mm_add_epi16(xmme128i, _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+8+ 2*AccumOffset]));
                            xmma128i=_mm_srli_si128(xmma128i, 4);
                            xmmb128i=_mm_slli_si128(xmme128i, 12);
                            xmmb128i=_mm_add_epi16(xmmb128i,xmma128i);
                            IQSum128i[Slice+2]=_mm_add_epi32( _mm_madd_epi16(xmmb128i,InterCoef[IndCoef  ]), IQSum128i[Slice+2]);
                            IQSum128i[Slice+3]=_mm_add_epi32( _mm_madd_epi16(xmmb128i,InterCoef[IndCoef+1]), IQSum128i[Slice+3]);
                            xmma128i=_mm_srli_si128(xmma128i, 4);
                            xmmc128i=_mm_slli_si128(xmme128i, 8);
                            xmmc128i=_mm_add_epi16(xmmc128i,xmma128i);
                            IQSum128i[Slice+4]=_mm_add_epi32( _mm_madd_epi16(xmmc128i,InterCoef[IndCoef  ]), IQSum128i[Slice+4]);
                            IQSum128i[Slice+5]=_mm_add_epi32( _mm_madd_epi16(xmmc128i,InterCoef[IndCoef+1]), IQSum128i[Slice+5]);
                            xmma128i=_mm_srli_si128(xmma128i, 4);
                            xmmd128i=_mm_slli_si128(xmme128i, 4);
                            xmmd128i=_mm_add_epi16(xmmd128i,xmma128i);
                            IQSum128i[Slice+6]=_mm_add_epi32( _mm_madd_epi16(xmmd128i,InterCoef[IndCoef  ]), IQSum128i[Slice+6]);
                            IQSum128i[Slice+7]=_mm_add_epi32( _mm_madd_epi16(xmmd128i,InterCoef[IndCoef+1]), IQSum128i[Slice+7]);
                            IQSum128i[Slice+8]=_mm_add_epi32( _mm_madd_epi16(xmme128i,InterCoef[IndCoef  ]), IQSum128i[Slice+8]);
                            IQSum128i[Slice+9]=_mm_add_epi32( _mm_madd_epi16(xmme128i,InterCoef[IndCoef+1]), IQSum128i[Slice+9]);
                           
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;

                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1610504 :       // 5 pixels computed per slice without folding at 100% BandWidth with accumulation of 4 firings and immediate lut
                // Start loop on channels
                for (k=0; k < NChannelsRight+NChannelsLeft; k++) {      // Loop on contributing piezos
                    if (k == NChannelsRight) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if (k1 == NChannelsRight+NChannelsLeft) k1 = k1 -1; // last Channel case
                    if (k1 < NChannelsRight) {
                        ChannelNext = Absc + k1;
                    }
                    else {
                        ChannelNext = Absc - (k1-NChannelsRight) -1;
                    }
                    
                    if (ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;    // Acquisition Channel
                    ChannelStartPtr= (short *) &RcvDataShort[(ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq))) +
                                                               4*NSamples*(1+SyntAcq)*NRec*(1-PlaneWave) +     // multiply by 4 because there are four firings per recon
                                                                NSamples*(Channel/(NChannels/(1+SyntAcq)))];
                    DeltaPrefetch = (
                                    (ChannelOffset/2)*(ChannelNext%(NChannels/(1+SyntAcq))) +       
                                    NSamples*(ChannelNext/(NChannels/(1+SyntAcq))) -
                                    (ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq))) -
                                    NSamples*(Channel/(NChannels/(1+SyntAcq)))
                                    )/8;
                    Channel = ChannelNext;
                    AccumOffset=NSamples*(1+SyntAcq);
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for (j=0; j<NLines; j++) {
                        SliceStart = RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for (Slice = (2*SliceStart + j*2*NSlices)*5; Slice <(2*NSlices + j*2*NSlices)*5; Slice += 10) {
                            Time32 = RecoLut[RecoLutInd++];          // get time value in lambda*32
                            IndCoef = RecoLut[RecoLutInd++];
                            dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32] + DeltaPrefetch);
                            _mm_prefetch(dataprefetch,  _MM_HINT_NTA);
                            dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32+AccumOffset] + DeltaPrefetch);
                            _mm_prefetch(dataprefetch,  _MM_HINT_NTA);
                            xmma128i =_mm_add_epi16(_mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 0]), _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ AccumOffset]));
                            xmma128i =_mm_add_epi16(xmma128i, _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 2*AccumOffset]));
                            xmma128i =_mm_add_epi16(xmma128i, _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 3*AccumOffset]));
                            IQSum128i[Slice]=_mm_add_epi32( _mm_madd_epi16(xmma128i,InterCoef[IndCoef  ]), IQSum128i[Slice+0]);
                            IQSum128i[Slice+1]=_mm_add_epi32( _mm_madd_epi16(xmma128i,InterCoef[IndCoef+1]), IQSum128i[Slice+1]);
                            xmme128i =_mm_add_epi16(_mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 8]), _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+8+ AccumOffset]));
                            xmme128i =_mm_add_epi16(xmme128i, _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+8+ 2*AccumOffset]));
                            xmme128i =_mm_add_epi16(xmme128i, _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+8+ 3*AccumOffset]));
                            xmma128i=_mm_srli_si128(xmma128i, 4);
                            xmmb128i=_mm_slli_si128(xmme128i, 12);
                            xmmb128i=_mm_add_epi16(xmmb128i,xmma128i);
                            IQSum128i[Slice+2]=_mm_add_epi32( _mm_madd_epi16(xmmb128i,InterCoef[IndCoef  ]), IQSum128i[Slice+2]);
                            IQSum128i[Slice+3]=_mm_add_epi32( _mm_madd_epi16(xmmb128i,InterCoef[IndCoef+1]), IQSum128i[Slice+3]);
                            xmma128i=_mm_srli_si128(xmma128i, 4);
                            xmmc128i=_mm_slli_si128(xmme128i, 8);
                            xmmc128i=_mm_add_epi16(xmmc128i,xmma128i);
                            IQSum128i[Slice+4]=_mm_add_epi32( _mm_madd_epi16(xmmc128i,InterCoef[IndCoef  ]), IQSum128i[Slice+4]);
                            IQSum128i[Slice+5]=_mm_add_epi32( _mm_madd_epi16(xmmc128i,InterCoef[IndCoef+1]), IQSum128i[Slice+5]);
                            xmma128i=_mm_srli_si128(xmma128i, 4);
                            xmmd128i=_mm_slli_si128(xmme128i, 4);
                            xmmd128i=_mm_add_epi16(xmmd128i,xmma128i);
                            IQSum128i[Slice+6]=_mm_add_epi32( _mm_madd_epi16(xmmd128i,InterCoef[IndCoef  ]), IQSum128i[Slice+6]);
                            IQSum128i[Slice+7]=_mm_add_epi32( _mm_madd_epi16(xmmd128i,InterCoef[IndCoef+1]), IQSum128i[Slice+7]);
                            IQSum128i[Slice+8]=_mm_add_epi32( _mm_madd_epi16(xmme128i,InterCoef[IndCoef  ]), IQSum128i[Slice+8]);
                            IQSum128i[Slice+9]=_mm_add_epi32( _mm_madd_epi16(xmme128i,InterCoef[IndCoef+1]), IQSum128i[Slice+9]);
                           
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x10504 :       // 5 pixels computed per slice without folding at 100% BandWidth
                #ifdef __ASMGCC64_OK
                if ( ReconMethod == E_ASM )
                {
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    __asm__ __volatile__
                    (
                        ////////////////////////////////////////////////////////////////////////////
                        // Loop on contributing piezos
                        " xor   %%r10,              %%r10           \n" //r10=0
                        "9:                                         \n" //main loop on channels
                        " xor   %%r14,              %%r14           \n" //r14=0
                        " xor   %%r15,              %%r15           \n" //r15=0
                        " mov   %%r10d,             %%r15d          \n" //get the channel counter
                        " cmp   %%r15d,             28(%[pASMLut])  \n" //if (k == NChannelsRight)
                        " cmove 40(%[pASMLut]),     %[RecoLutInd]   \n" //RecoLutInd = RecoLutMoinsStart[0]

                        " inc   %%r15d                              \n" //k1 = k +1
                        " cmp   %%r15d,             32(%[pASMLut])  \n" //if (k1 == NChannelsRight+NChannelsLeft)
                        " jne   11f                                 \n"
                        " dec   %%r15d                              \n" //k1 = k1 -1
                        "11:                                        \n"
                        " mov   48(%[pASMLut]),     %%r14d          \n"
                        " cmp   28(%[pASMLut]),     %%r15d          \n" //if (k1 < NChannelsRight)
                        " jb    12f                                 \n" //jump to 12
                        " sub   %%r15d,             %%r14d          \n"
                        " sub   %%r15d,             %%r14d          \n"
                        " add   28(%[pASMLut]),     %%r14d          \n"
                        " dec   %%r14d                              \n"
                        "12:                                        \n"
                        " add   %%r15d,             %%r14d          \n"

                        " cmp   36(%[pASMLut]),     %%r14d          \n" //if ChannelNext < NChannel,
                        " jb    10f                                 \n" //jump after the next instruction
                        " sub   36(%[pASMLut]),     %%r14d          \n" // else ChannelNext = ChannelNext - NChannels
                        "10:                                        \n"
                        " xor   %%r13,          %%r13  \n" //r13=0
                        " movl  %[Channel],     %%eax  \n" //get Channel ==> low dividend part
                        " xor   %%edx,          %%edx  \n" //high dividend part is null
                        " divl 12(%[pASMLut])          \n" //channel % get pASMLut[3] => remainder in edx
                        " mov   %%edx,          %%eax  \n"
                        " mull  8(%[pASMLut])          \n" //multiply result by pASMLut[2] => result in eax
                        " mov  16(%[pASMLut]),  %%r13d \n" //get pASMLut[4]
                        " add   %%eax,          %%r13d \n" //add pASMLut[4]+previous result ==> in r13

                        " xor   %%edx,          %%edx  \n" //high dividend part is null
                        " movl  %[Channel],     %%eax  \n" //get Channel ==> low dividend part
                        " divl  12(%[pASMLut])         \n" //channel  get pASMLut[5] => quotient in eax
                        " mull  24(%[pASMLut])         \n" //multiply result by pASMLut[6] (NSamples) => result in eax
                        " add   %%eax,          %%r13d \n" //both previous result ==> in r13

                        " mov   %[RcvDataShort], %[ChannelStartPtr]\n" //and create the final RcvDataShortPtr
                        " add   %%r13, %[ChannelStartPtr]\n"
                        " add   %%r13, %[ChannelStartPtr]\n" //twice, because it's shorts

                        " mov   %%r14d,         %[Channel] \n"  //Channel = ChannelNext

                        " xor   %%r13,      %%r13 \n" //r13=0 //the lines counter
//                         " prefetch 64(%[RecoLut], %%rbx, 2) \n" //prefetch the lut to be keep on L2 after use

                        ////////////////////////////////////////////////////////////////////////////
                        "3:\n " // line loop
                        " xor   %%r14,      %%r14 \n" //r14=0
                        " xor   %%r8,      %%r8   \n" //r8=0

                        //compute index values
                        //compute slice start for r14
                        " mov  (%[RecoLut], %%rbx, 2),     %%r14w \n" //get SliceStart (to replace after by SliceStart = RecoLut[RecoLutInd++];

                        " lea  (%%r14, %%r14),     %%r14 \n" //= SliceStart*2
                        " mov  (%[pASMLut]),       %%r8d \n" //get NSlices
                        " lea  (%%r8, %%r8),       %%r8  \n" //= NSlices*2
                        " mov   %%r13,             %%r9  \n" //get j
                        " imul  %%r8,              %%r9  \n" //= j*2*NSlices
                        " add   %%r9,              %%r14 \n" //= SliceStart*2 + j*2*NSlices
                        " add   %%r8,              %%r9  \n" //= j*2*NSlices + NSlices*2
                        " mov   %%r14,             %%r8  \n" //get a slice_start copy
                        " sub   %%r8,              %%r9  \n" //= slice_start - slice_stop ==> nb_loop*2
                        " lea  (%%r14, %%r14,4),   %%r14 \n" //= (SliceStart*2 + j*2*NSlices)*5    ==> slice_start
                        " lea  (%%r14, %%r14),     %%r14 \n" //= (SliceStart*2 + j*2*NSlices)*5*2
                        " lea  (      , %%r14,8),  %%r14 \n" //= (SliceStart*2 + j*2*NSlices)*5*2*8 ==> offset for &IQSum128i[slice_start];
                        //" lea  (%%r9d, %%r9d,2),   %%r9d  \n" //= (j*2*NSlices + 2*NSlices)*3 ==> slice_stop

                        //validate values and start loop
                        " mov %%r9,                %%r15  \n"
                        " add %[RecoLutInd],       %%r15d \n"

                        " inc                 %[RecoLutInd]   \n" //nb_loop_counter ++

                        " cmp   $0,                %%r9   \n" //compare loop value with 0, to do not enter in loop if not necessary
                        " jle 2f                          \n" //immediat output case

                        ////////////////////////////////////////////////////////////////////////////
                        "1:\n"  //slice loop
                        //prefetch the lut
                        //load and compute corrects IndCoef and Time32
                        " mov    %[Delay32],            %%r9d  \n" //Load Delay32
                        " addw  (%[RecoLut], %%rbx, 2), %%r9w  \n" //Time 32 = Delay32 + Time32

                        " mov       $0x1E,              %%rax  \n" //load 0x1E for IndCoef
                        " and       %%r9,               %%rax  \n" //IndCoef = (Time32)& 0x1E

                        " cmp       $7,                 %%rax  \n" //if IndCoef < $7
                        " setle                         %%dl   \n" // prepare the value in rdx to decrease Time32

                        " inc    %[RecoLutInd]                 \n" //nb_loop_counter ++
                        " add   (%[RecoLut], %%rbx, 2),  %%ax  \n" //Load Lut IndCoef
                        " shl     $4,                    %%ax  \n" //InterCoef = InterCoef *4 Because InterCoef type is xmm128

//                         " prefetch  16(%[InterCoef], %%rax) \n"


                        //shift Time 32 and mask it
                        " shr    $4,               %%r9    \n" //Time32 = Tim3e2 >> 4
                        " and    $0xFFE,           %%r9    \n" //Time32 = Time32 & 0xFFE

                        " movdqa (%[InterCoef], %%rax),      %%xmm8  \n" //load InterCoef from    pRecolut + IndCoef

                        //and decrease it following the previous result in rdx
                        " sub    %%rdx,             %%r9   \n" //if (Indcoef > 7) Time32 --

                        //get data[0]
                        " movdqu (%[ChannelStartPtr], %%r9, 2),    %%xmm0  \n" //load the source value from pSrc + Time32*sizeof(short)

                        //Compute slice[0]:
                        " movdqa %%xmm8,                    %%xmm7  \n" //use positions dcba
                        " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between source and lut
                        " paddd  (%[IQSum128i], %%r14),     %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm7,     (%[IQSum128i],%%r14)  \n" //copy result to dest

                            " movdqa 16(%[InterCoef], %%rax),    %%xmm9  \n" //load the InterCoef from pRecolut + IndCoef+16 (to get the next 128 block)

                        //Compute slice[1]:
                        " movdqa %%xmm9,                    %%xmm7  \n"
                        " pmaddwd %%xmm0,                   %%xmm7  \n" //_mm_madd_epi16 between source and lut
                        " paddd  16(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest[Slice+1]
                        " movdqa  %%xmm7,   16(%[IQSum128i],%%r14)  \n" //copy result to Dest[Slice+1]

                            //get data[1]
                            " movdqu 16(%[ChannelStartPtr], %%r9, 2),  %%xmm1  \n" //load the source value from pSrc + Time32*sizeof(short) + 16

                        //Compute slice[8]:
                        " movdqa  %%xmm8,                   %%xmm7  \n" // use positions hgfe
                        " pmaddwd %%xmm1,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  128(%[IQSum128i], %%r14),  %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm7,  128(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //Compute slice[9]:
                        " movdqa  %%xmm9,                   %%xmm7  \n"
                        " pmaddwd %%xmm1,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  144(%[IQSum128i], %%r14),  %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm7,  144(%[IQSum128i],%%r14)  \n" //copy result to dest*/

                        //intermediate calculs
                        " pshufd $0x4E, %%xmm0,             %%xmm2  \n"
                        " punpcklqdq    %%xmm1,             %%xmm2  \n" //create the position for slides 4 & 5 => fedc

                        //Compute slice[4]:
                        " movdqa %%xmm8,                    %%xmm7  \n"
                        " pmaddwd %%xmm2,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  64(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm7,   64(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //Compute slice[5]:
                        " pmaddwd %%xmm9,                   %%xmm2  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  80(%[IQSum128i], %%r14),   %%xmm2  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm2,   80(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //intermediate calculs
                        " movss %%xmm1,                     %%xmm0 \n"  //create the position for slices 2 && 3 => edcb
                        " pshufd $0x39, %%xmm0,             %%xmm2 \n"

                        //Compute slice[2]:
                        " movdqa %%xmm8,                    %%xmm7  \n"
                        " pmaddwd %%xmm2,                   %%xmm7  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  32(%[IQSum128i], %%r14),   %%xmm7  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm7,   32(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //Compute slice[3]:
                        " pmaddwd %%xmm9,                   %%xmm2  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  48(%[IQSum128i], %%r14),   %%xmm2  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm2,   48(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //intermediate calculs
                        " pshufd $0xFF, %%xmm0,             %%xmm0 \n"  //move the old MSB to LSB
                        " pslldq  $4,                       %%xmm1 \n"
                        " movss %%xmm0,                     %%xmm1 \n"  //create the position for slices 6 && 7 => gfed

                        //Compute slice[6]:
                        " pmaddwd %%xmm1,                   %%xmm8  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  96(%[IQSum128i], %%r14),   %%xmm8  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm8,   96(%[IQSum128i],%%r14)  \n" //copy result to dest

                        //Compute slice[7]:
                        " pmaddwd %%xmm9,                   %%xmm1  \n" //_mm_madd_epi16 between result of intermediate part and lut
                        " paddd  112(%[IQSum128i], %%r14),  %%xmm1  \n" //_mm_add_epi32 between previous result and Dest
                        " movdqa  %%xmm1,  112(%[IQSum128i],%%r14)  \n" //copy result to dest

                        " inc                        %[RecoLutInd]  \n" //nb_loop_counter ++

                        " add     $160,                      %%r14   \n" //increase slice pointer (i128*10 = 160)

                        " cmp     %%r15d,            %[RecoLutInd]  \n" //nb_loop_counter = nb_loop
                        " jb       1b                               \n" //and branche it if non zero equal
                        "2:                                             \n" //on se casse...
                        " inc %%r13 \n " //increment line loop
                        " cmp  4(%[pASMLut]), %%r13d \n " //if < NLines
                        " jb 3b \n " //jump to loop beginning

                        " inc                       %%r10d          \n"
                        " cmp       %%r10d,         32(%[pASMLut])  \n"
                        " jne       9b                              \n"//*/

                    :   [IQSum128i]         "+r" ( IQSum128i ),
                        [ChannelStartPtr]   "+r" ( ChannelStartPtr ),
                        [RecoLutInd]        "+b" ( RecoLutInd ),
                        [RecoLut]           "+r" ( RecoLut ),
                        [InterCoef]         "+r" ( InterCoef ),
                        [pASMLut]           "+r" ( pASMLut ),
                        [Channel]           "+m" ( Channel ),
                        [ChannelNext]       "+m" ( ChannelNext )
                        :   [RcvDataShort]      "m" ( RcvDataShort ),
                        [Delay32]           "m" ( Delay32 )

                    :   "%r15", "%r14", "%r13", "%r12", "%r10", "%r9", "%r8", "%rax", "%rdx",
                        "%xmm0", "%xmm1", "%xmm2", "%xmm3", "%xmm7", "%xmm8", "%xmm9"
                    );
                } //end ASM
                else
                #endif //end #ifdef __ASMGCC64_OK
                {
                    // Start loop on channels
                    for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                    {
                        if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                        // Compute next Channel for prefetch
                        k1 = k + 1;
                        if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                        if ( k1 < NChannelsRight )
                        {
                            ChannelNext = Absc + k1;
                        }
                        else
                        {
                            ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                        }

                        if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                        ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                        DeltaPrefetch = (
                                            ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                            NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                            ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                            NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                        ) /8;

                        Channel = ChannelNext;
                        ////////////////////////////////////////////////////////////////////////////////
                        // Loop on lines
                        for ( j=0; j<NLines; j++ )
                        {
                            SliceStart = RecoLut[RecoLutInd++];
                            ////////////////////////////////////////////////////////////////////////////
                            //loop on the Slices for this region
                            for ( Slice = ( 2*SliceStart + j*2*NSlices ) *5; Slice < ( 2*NSlices + j*2*NSlices ) *5; Slice += 10 )
                            {
                                Time32 = RecoLut[RecoLutInd++];          // get time value in lambda*32
                                Time32 = Time32 + Delay32;               // Time32 corrected by delay
                                IndCoef = ( Time32 & 0x001E );//+ RecoLut[RecoLutInd++];
                                if ( IndCoef >7 )           // Phase is >=4
                                {
                                    Time32 = ( Time32>>4 ) & 0xFFE;
                                    IndCoef = IndCoef+ RecoLut[RecoLutInd++];
                                    dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                                    _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                                    xmma128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] );
                                    IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                                    IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                                    xmme128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] );
                                    xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                    xmmb128i=_mm_slli_si128 ( xmme128i, 12 );
                                    xmmb128i=_mm_add_epi16 ( xmmb128i,xmma128i );
                                    IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                                    IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                                    xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                    xmmc128i=_mm_slli_si128 ( xmme128i, 8 );
                                    xmmc128i=_mm_add_epi16 ( xmmc128i,xmma128i );
                                    IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                                    IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                                    xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                    xmmd128i=_mm_slli_si128 ( xmme128i, 4 );
                                    xmmd128i=_mm_add_epi16 ( xmmd128i,xmma128i );
                                    IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                                    IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                                    IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                                    IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );
                                }
                                else                        // Phase is between 0 and 3
                                {
                                    Time32 = ( ( Time32>>4 ) & 0xFFE ) -1;
                                    IndCoef = IndCoef+ RecoLut[RecoLutInd++];
                                    dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                                    _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                                    xmma128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] );
                                    IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                                    IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                                    xmme128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] );
                                    xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                    xmmb128i=_mm_slli_si128 ( xmme128i, 12 );
                                    xmmb128i=_mm_add_epi16 ( xmmb128i,xmma128i );
                                    IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                                    IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                                    xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                    xmmc128i=_mm_slli_si128 ( xmme128i, 8 );
                                    xmmc128i=_mm_add_epi16 ( xmmc128i,xmma128i );
                                    IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                                    IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                                    xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                    xmmd128i=_mm_slli_si128 ( xmme128i, 4 );
                                    xmmd128i=_mm_add_epi16 ( xmmd128i,xmma128i );
                                    IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                                    IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                                    IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                                    IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );
                                }
                            }   // End loop on slices
                        }       // End loop on Lines
                    }           // End loop on Channels
                } //end Intrinsics
                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x0210504 :       // 5 pixels computed per slice without folding at 100% BandWidth with accumulation
                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }

                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     2*NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +      //multiply by two because there are two firings per recon
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    AccumOffset=NSamples* ( 1+SyntAcq );
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart = RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ) *5; Slice < ( 2*NSlices + j*2*NSlices ) *5; Slice += 10 )
                        {
                            Time32 = RecoLut[RecoLutInd++];          // get time value in lambda*32
                            Time32 = Time32 + Delay32;               // Time32 corrected by delay
                            IndCoef = ( Time32 & 0x001E );//+ RecoLut[RecoLutInd++];
                            if ( IndCoef >7 )           // Phase is >=4
                            {
                                Time32 = ( Time32>>4 ) & 0xFFE;
                                IndCoef = IndCoef+ RecoLut[RecoLutInd++];
                                dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                                _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                                dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32+AccumOffset] + DeltaPrefetch );
                                _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                                //xmma128i = _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 0]);
                                xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + AccumOffset] ) );
                                IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                                IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                                //xmme128i = _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 8]);
                                xmme128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) +8+ AccumOffset] ) );
                                xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                xmmb128i=_mm_slli_si128 ( xmme128i, 12 );
                                xmmb128i=_mm_add_epi16 ( xmmb128i,xmma128i );
                                IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                                IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                                xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                xmmc128i=_mm_slli_si128 ( xmme128i, 8 );
                                xmmc128i=_mm_add_epi16 ( xmmc128i,xmma128i );
                                IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                                IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                                xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                xmmd128i=_mm_slli_si128 ( xmme128i, 4 );
                                xmmd128i=_mm_add_epi16 ( xmmd128i,xmma128i );
                                IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                                IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                                IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                                IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );
                            }
                            else                        // Phase is between 0 and 3
                            {
                                Time32 = ( ( Time32>>4 ) & 0xFFE ) -1;
                                IndCoef = IndCoef+ RecoLut[RecoLutInd++];
                                //dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32] + DeltaPrefetch);
                                //_mm_prefetch(dataprefetch,  _MM_HINT_NTA);
                                //xmma128i = _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 0]);
                                xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + AccumOffset] ) );
                                IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                                IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                                //xmme128i = _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32)+ 8]);
                                xmme128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) +8+ AccumOffset] ) );
                                xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                xmmb128i=_mm_slli_si128 ( xmme128i, 12 );
                                xmmb128i=_mm_add_epi16 ( xmmb128i,xmma128i );
                                IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                                IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                                xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                xmmc128i=_mm_slli_si128 ( xmme128i, 8 );
                                xmmc128i=_mm_add_epi16 ( xmmc128i,xmma128i );
                                IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                                IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                                xmma128i=_mm_srli_si128 ( xmma128i, 4 );
                                xmmd128i=_mm_slli_si128 ( xmme128i, 4 );
                                xmmd128i=_mm_add_epi16 ( xmmd128i,xmma128i );
                                IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                                IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                                IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                                IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );
                            }
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;

                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x10100 :       // Only one pixel computed per slice without folding at 100% BandWidth
                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }

                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    /*    DeltaPrefetch = (
                    (ChannelOffset/2)*(ChannelNext%(NChannels/(1+SyntAcq))) +
                    NSamples*(ChannelNext/(NChannels/(1+SyntAcq))) -
                    (ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq))) -
                    NSamples*(Channel/(NChannels/(1+SyntAcq)))
                    )/8; */
                    Channel = ChannelNext;
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart = RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ); Slice < ( 2*NSlices + j*2*NSlices ); Slice += 2 )
                        {
                            Time32 = RecoLut[RecoLutInd++];          // get time value in lambda*32
                            Time32 = Time32 + Delay32;               // Time32 corrected by delay
                            IndCoef = ( Time32 & 0x001E );//+ RecoLut[RecoLutInd++];
                            if( RecoLutInd >= recon_lut_size_nb )
                            {
                              printf("Error: RecoLutInd >= recon_lut_size_nb\n");
                              printf("recon_lut_size_nb = %d, RecoLutInd = %d\n", recon_lut_size_nb, RecoLutInd);
                              return 0;
                            }

                            if ( IndCoef >7 )           // Phase is >=4
                            {
                                Time32 = ( Time32>>4 ) & 0xFFE;
                                IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[Time32] ),InterCoef[IndCoef + RecoLut[RecoLutInd] ] ), IQSum128i[Slice  ] );
                                IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[Time32] ),InterCoef[IndCoef+ RecoLut[RecoLutInd++]+1] ), IQSum128i[Slice+1] );
                            }
                            else                        // Phase is between 0 and 3
                            {
                                Time32 = ( ( Time32>>4 ) & 0xFFE ) -1;
                                IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[Time32] ),InterCoef[IndCoef + RecoLut[RecoLutInd] ] ), IQSum128i[Slice  ] );
                                IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[Time32] ),InterCoef[IndCoef+ RecoLut[RecoLutInd++]+1] ), IQSum128i[Slice+1] );
                            }
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels


                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1120504 :       // 5 pixels computed per slice with folding at 200% BandWidth and immediate Lut. Same code for 0x1110508
                // Determine number of folded channels
                if ( ( NChannelsLeft+CentralChannel ) <= NChannelsRight )
                {
                    MinRightLeft = NChannelsLeft+CentralChannel;           // Number of folded channels
                    MaxRightLeft = NChannelsRight;          // Max Channel number on one side
                }
                else
                {
                    MaxRightLeft = NChannelsLeft+CentralChannel;
                    MinRightLeft = NChannelsRight;
                }
                ////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////
                // Case where there is a central channel
                if ( CentralChannel ==1 )
                {
                    Channel= Absc;
                    // Acquisition Channel
                    if ( Channel >= NChannels ) Channel = Channel - NChannels;   // Acquisition Channel right
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    /*                     DeltaPrefetch = (
                    (ChannelOffset/2)*(ChannelNext%(NChannels/(1+SyntAcq))) +
                    NSamples*(ChannelNext/(NChannels/(1+SyntAcq))) -
                    (ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq))) -
                    NSamples*(Channel/(NChannels/(1+SyntAcq)))
                    )/8;*/
                    ////////////////////////////////////////////////////////////////////////////////
                    // Start the line
                    SliceStart =RecoLut[RecoLutInd++];
                    //loop on the Slices for this region
                    for ( Slice = ( 2*SliceStart ) *5; Slice < ( 2*NSlices ) *5; Slice += 10 )
                    {
                        Time32 =RecoLut[RecoLutInd++];
                        IndCoef =RecoLut[RecoLutInd++];
                        //dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32] + DeltaPrefetch);
                        //_mm_prefetch(dataprefetch,  _MM_HINT_NTA);
                        xmma128i =_mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] );
                        IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                        IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                        xmmc128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] );
                        xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                        xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                        IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                        IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                        IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                        IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                        xmme128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 16] );
                        xmmd128i = _mm_unpacklo_epi64 ( xmmd128i, xmme128i );
                        IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                        IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                        IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                        IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );
                    }   // End loop on slices
                }           // End loop on Channels

                ////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////
                // Loop on folded channels
                for ( k=CentralChannel; k < MinRightLeft; k++ )
                {
                    ChannelRight = Absc + k;
                    ChannelLeft = Absc - k -1+ CentralChannel;
                    // Acquisition Channel
                    if ( ChannelRight >= NChannels ) ChannelRight = ChannelRight - NChannels;   // Acquisition Channel right
                    if ( ChannelLeft >= NChannels ) ChannelLeft = ChannelLeft - NChannels;   // Acquisition Channel left
                    ChannelStartPtrRight= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( ChannelRight% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                          NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                          NSamples* ( ChannelRight/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    ChannelStartPtrLeft= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( ChannelLeft% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                         NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                         NSamples* ( ChannelLeft/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    /*                     DeltaPrefetch = (
                    (ChannelOffset/2)*(ChannelNext%(NChannels/(1+SyntAcq))) +
                    NSamples*(ChannelNext/(NChannels/(1+SyntAcq))) -
                    (ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq))) -
                    NSamples*(Channel/(NChannels/(1+SyntAcq)))
                    )/8; */
                    ////////////////////////////////////////////////////////////////////////////////
                    // Start line
                    SliceStart =RecoLut[RecoLutInd++];
                    //loop on the Slices for this region
                    for ( Slice = ( 2*SliceStart ) *5; Slice < ( 2*NSlices ) *5; Slice += 10 )
                    {
                        Time32 =RecoLut[RecoLutInd++];
                        IndCoef =RecoLut[RecoLutInd++];
                        //dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32] + DeltaPrefetch);
                        //_mm_prefetch(dataprefetch,  _MM_HINT_NTA);
                        xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtrRight[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtrLeft[ ( Time32 ) + 0] ) );
                        IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                        IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                        xmmc128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtrRight[ ( Time32 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtrLeft[ ( Time32 ) + 8] ) );
                        xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                        xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                        IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                        IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                        IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                        IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                        xmme128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtrRight[ ( Time32 ) +16] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtrLeft[ ( Time32 ) +16] ) );
                        xmmd128i = _mm_unpacklo_epi64 ( xmmd128i, xmme128i );
                        IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                        IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                        IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                        IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );
                    }   // End loop on slices
                }           // End loop on Channels
                ////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////
                // Loop on unfolded channels
                for ( k=MinRightLeft; k < MaxRightLeft; k++ )
                {
                    if ( NChannelsRight > NChannelsLeft )
                    {
                        Channel = Absc + k;
                    }
                    else
                    {
                        Channel= Absc - k -1+ CentralChannel;
                    }
                    // Acquisition Channel
                    if ( Channel >= NChannels ) Channel = Channel - NChannels;   // Acquisition Channel right
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    /*                     DeltaPrefetch = (
                    (ChannelOffset/2)*(ChannelNext%(NChannels/(1+SyntAcq))) +
                    NSamples*(ChannelNext/(NChannels/(1+SyntAcq))) -
                    (ChannelOffset/2)*(Channel%(NChannels/(1+SyntAcq))) -
                    NSamples*(Channel/(NChannels/(1+SyntAcq)))
                    )/8; */
                    ////////////////////////////////////////////////////////////////////////////////
                    // start line
                    SliceStart = RecoLut[RecoLutInd++];
                    //loop on the Slices for this region
                    for ( Slice = ( 2*SliceStart ) *5; Slice < ( 2*NSlices ) *5; Slice += 10 )
                    {
                        Time32 =RecoLut[RecoLutInd++];
                        IndCoef =RecoLut[RecoLutInd++];
                        //dataprefetch= (char *) ((__m128i *) &ChannelStartPtr[Time32] + DeltaPrefetch);
                        //_mm_prefetch(dataprefetch,  _MM_HINT_NTA);
                        xmma128i =_mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] );
                        IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                        IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                        xmmc128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] );
                        xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                        xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                        IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                        IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                        IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                        IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                        xmme128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 16] );
                        xmmd128i = _mm_unpacklo_epi64 ( xmmd128i, xmme128i );
                        IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                        IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                        IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                        IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );
                    }   // End loop on slices
                }           // End loop on Channels
                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1220504 :       // 5 pixels computed per slice without folding at 200% BandWidth with accumulation and immediate Lut


                xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + AccumOffset] ) );


                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }

                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     2*NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +     // multiply by 2 because there are two firings per recon
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    AccumOffset=NSamples* ( 1+SyntAcq );
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart =RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ) *5; Slice < ( 2*NSlices + j*2*NSlices ) *5; Slice += 10 )
                        {
                            Time32 = RecoLut[RecoLutInd++];
                            IndCoef=RecoLut[RecoLutInd++];
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32+AccumOffset] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + AccumOffset] ) );
                            IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            xmmc128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8 +AccumOffset] ) );
                            xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                            xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                            xmme128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 16] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 16 + AccumOffset] ) );
                            xmmd128i = _mm_unpacklo_epi64 ( xmmd128i, xmme128i );
                            IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                            IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                            IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                            IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;
                                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1420504 :       // 5 pixels computed per slice without folding at 200% BandWidth with accumulation of 3 firings and immediate Lut


                // xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + AccumOffset] ) );


                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }

                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     3*NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +     // multiply by 2 because there are two firings per recon
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    AccumOffset=NSamples* ( 1+SyntAcq );
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart =RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ) *5; Slice < ( 2*NSlices + j*2*NSlices ) *5; Slice += 10 )
                        {
                            Time32 = RecoLut[RecoLutInd++];
                            IndCoef=RecoLut[RecoLutInd++];
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32+AccumOffset] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + AccumOffset] ) );
                            xmma128i =_mm_add_epi16 (xmma128i, _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 2*AccumOffset] ) );
                            IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            xmmc128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8 +AccumOffset] ) );
                            xmmc128i =_mm_add_epi16 ( xmmc128i, _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8 +2*AccumOffset] ) );
                            xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                            xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                            xmme128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 16] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 16 + AccumOffset] ) );
                            xmme128i =_mm_add_epi16 (xmme128i, _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 16 + 2*AccumOffset] ) );
                            xmmd128i = _mm_unpacklo_epi64 ( xmmd128i, xmme128i );
                            IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                            IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                            IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                            IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x1020504 :       // 5 pixels computed per slice without folding at 200% BandWidth and immediate Lut
                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }

                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart =RecoLut[RecoLutInd++];

                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ) *5; Slice < ( 2*NSlices + j*2*NSlices ) *5; Slice += 10 )
                        {
                            Time32 = RecoLut[RecoLutInd++];
                            IndCoef=RecoLut[RecoLutInd++];
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            xmma128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 0] );
                            IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            xmmc128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 8] );
                            xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                            xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                            xmme128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 ) + 16] );
                            xmmd128i = _mm_unpacklo_epi64 ( xmmd128i, xmme128i );
                            IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                            IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                            IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                            IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x20504 :       // Five pixels computed at a time per slice, with lambda/2 pitch
                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }
                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart = RecoLut[RecoLutInd++];
                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ) *5; Slice < ( 2*NSlices + j*2*NSlices ) *5; Slice += 10 )
                        {
                            Time32 = RecoLut[RecoLutInd++];          // get time value in lambda*32
                            Time32 = Time32 + Delay32;               // Time32 corrected by delay
                            IndCoef = ( Time32 & 0x0006 ) + RecoLut[RecoLutInd++];
                            //data = (__m128i *) &ChannelStartPtr[Time32 >> 3];
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32 >> 3] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            xmma128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 >> 3 ) + 0] );
                            IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            xmmc128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 >> 3 ) + 8] );
                            xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                            xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                            xmme128i = _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32 >> 3 ) + 16] );
                            xmmd128i = _mm_unpacklo_epi64 ( xmmd128i, xmme128i );
                            IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                            IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                            IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                            IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x0220504 :       // Five pixels computed at a time per slice, with lambda/2 pitch and with accumulation
                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }
                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     2*NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +     // multiply by 2 because there are 2 firings per recon
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    AccumOffset=NSamples* ( 1+SyntAcq );
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart = RecoLut[RecoLutInd++];
                        ////////////////////////////////////////////////////////////////////////////
                        //loop on the Slices for this region
                        for ( Slice = ( 2*SliceStart + j*2*NSlices ) *5; Slice < ( 2*NSlices + j*2*NSlices ) *5; Slice += 10 )
                        {
                            Time32 = RecoLut[RecoLutInd++];          // get time value in lambda*32
                            Time32 = Time32 + Delay32;               // Time32 corrected by delay
                            IndCoef = ( Time32 & 0x0006 ) + RecoLut[RecoLutInd++];
                            //data = (__m128i *) &ChannelStartPtr[Time32 >> 3];
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[Time32 >> 3] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            dataprefetch= ( char * ) ( ( __m128i * ) &ChannelStartPtr[ ( Time32>>3 ) +AccumOffset] + DeltaPrefetch );
                            _mm_prefetch ( dataprefetch,  _MM_HINT_NTA );
                            xmma128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32>>3 ) + 0] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32>>3 ) + AccumOffset] ) );
                            //xmma128i = _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32 >> 3)+ 0]);
                            IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                            IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( xmma128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            //xmmc128i = _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32 >> 3)+ 8]);
                            xmmc128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32>>3 ) + 8] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32>>3 ) + 8 +AccumOffset] ) );
                            xmmd128i = _mm_shuffle_epi32 ( xmmc128i, _MM_SHUFFLE ( 1,0,3,2 ) );
                            xmmb128i = _mm_unpackhi_epi64 ( xmma128i, xmmd128i );
                            IQSum128i[Slice+2]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+2] );
                            IQSum128i[Slice+3]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmb128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+3] );
                            IQSum128i[Slice+4]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+4] );
                            IQSum128i[Slice+5]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmc128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+5] );
                            //xmme128i = _mm_loadu_si128((__m128i *) &ChannelStartPtr[(Time32 >> 3)+ 16]);
                            xmme128i =_mm_add_epi16 ( _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32>>3 ) + 16] ), _mm_loadu_si128 ( ( __m128i * ) &ChannelStartPtr[ ( Time32>>3 ) + 16 + AccumOffset] ) );
                            xmmd128i = _mm_unpacklo_epi64 ( xmmd128i, xmme128i );
                            IQSum128i[Slice+6]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+6] );
                            IQSum128i[Slice+7]=_mm_add_epi32 ( _mm_madd_epi16 ( xmmd128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+7] );
                            IQSum128i[Slice+8]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef  ] ), IQSum128i[Slice+8] );
                            IQSum128i[Slice+9]=_mm_add_epi32 ( _mm_madd_epi16 ( xmme128i,InterCoef[IndCoef+1] ), IQSum128i[Slice+9] );
                        }   // End loop on slices
                    }       // End loop on Lines
                }           // End loop on Channels
                break;

                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////////////////////
            case  0x20100 :       // Only one pixel computed per slice without folding
                // Start loop on channels
                for ( k=0; k < NChannelsRight+NChannelsLeft; k++ )      // Loop on contributing piezos
                {
                    if ( k == NChannelsRight ) RecoLutInd = RecoLutMoinsStart[0]; // start reconstructing left
                    // Compute next Channel for prefetch
                    k1 = k + 1;
                    if ( k1 == NChannelsRight+NChannelsLeft ) k1 = k1 -1; // last Channel case
                    if ( k1 < NChannelsRight )
                    {
                        ChannelNext = Absc + k1;
                    }
                    else
                    {
                        ChannelNext = Absc - ( k1-NChannelsRight ) -1;
                    }
                    if ( ChannelNext >= NChannels ) ChannelNext = ChannelNext - NChannels;   // Acquisition Channel
                    ChannelStartPtr= ( short * ) &RcvDataShort[ ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                     NSamples* ( 1+SyntAcq ) *NRec* ( 1-PlaneWave ) +
                                     NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) ) ];
                    DeltaPrefetch = (
                                        ( ChannelOffset/2 ) * ( ChannelNext% ( NChannels/ ( 1+SyntAcq ) ) ) +
                                        NSamples* ( ChannelNext/ ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        ( ChannelOffset/2 ) * ( Channel% ( NChannels/ ( 1+SyntAcq ) ) ) -
                                        NSamples* ( Channel/ ( NChannels/ ( 1+SyntAcq ) ) )
                                    ) /8;
                    Channel = ChannelNext;
                    ////////////////////////////////////////////////////////////////////////////////
                    // Loop on lines
                    for ( j=0; j<NLines; j++ )
                    {
                        SliceStart = RecoLut[RecoLutInd++];
                        #ifdef __ASMGCC64_OK
                        if ( ReconMethod == E_ASM )
                        {
                            int slice_start = ( 2*SliceStart + j*2*NSlices );
                            int slice_stop  = ( 2*NSlices + j*2*NSlices )-1;
                            int nb_loop     = slice_stop - slice_start;
                            assert( slice_stop - slice_start > 0 );

                            short *pRecoLut = &RecoLut[RecoLutInd];
                            short *pSrc     = ChannelStartPtr;
                            __m128i *pDest  = &IQSum128i[slice_start];
                            int *pDestDebug = ( int * ) pDest;

                            long long _addr = 0; //for debug

                            __asm__ __volatile__
                            (
                                //initialisation
                                " xorq     %%r15,         %%r15    \n" //r15=0
                                " movl      %3,           %%r15d   \n" //counter for nb_loop
                                " mov      %%r15,         %%r14    \n" //get counter for Slices after
                                " dec      %%r14                   \n" //decrease loop
                                " shl       $4,           %%r14    \n" //and mult. it by 0x10 to create indirection for pDest (128bits/8)

                                "loop_20100:\n"

                                " xor     %%r8,            %%r8     \n" //r8 = 0
                                " mov    (%1, %%r15, 2),   %%r8w   \n" //IndCoef

                                " dec    %%r15d                     \n" //decrease loop
                                " movl    %7,              %%r9d    \n" //Delay32
                                " addw  (%1, %%r15, 2),    %%r9w    \n" //Delay32 += Time32

                                " mov    $0x6,             %%r10   \n"
                                " and   %%r9,              %%r10   \n" //(Delay32+Time32)& 0x06
                                " add   %%r10,             %%r8    \n" //IndCoef += (Dealy32+Time32)& 0x06

                                " shl    $4,               %%r8     \n" //Because InterCoef type is xmm128

                                " shr    $3,               %%r9    \n" //Time32 = (Time32+Delay32) >> 3

                                //Compute slice 0:
                                " movdqa (%4, %%r8),       %%xmm0   \n" //load the InterCoef value from pSrc + IndCoef*sizeof(short)
                                " movdqu (%5, %%r9, 2),    %%xmm1   \n" //load the source value from pSrc + Time32*sizeof(short)


                                " pmaddwd %%xmm1,          %%xmm0   \n" //_mm_madd_epi16 between source and lut
                                " paddd  (%0, %%r14),      %%xmm0   \n" //_mm_add_epi32 between previous result and Dest
                                " movdqa  %%xmm0,     (%0, %%r14)   \n" //copy result to dest

                                //Compute slice 1:
                                " movdqa 16(%4, %%r8),     %%xmm0  \n" //load the InterCoef value from pSrc + (IndCoef+1)*sizeof(short)
                                " pmaddwd %%xmm1,          %%xmm0  \n" //_mm_madd_epi16 between source and lut
                                " paddd  16(%0, %%r14),    %%xmm0  \n" //_mm_add_epi32 between previous result and Dest[Slice+1]
                                " movdqa  %%xmm0,   16(%0, %%r14)  \n" //copy result to Dest[Slice+1]

                                //It's finished, decrement counter and loop
                                " sub   $32,               %%r14  \n" //decrease loop
                                " dec                      %%r15  \n" //decrease loop
                                " jns    loop_20100               \n" //and branche it if non zero equal
                            : "+r" ( pDestDebug ), "+r" ( pRecoLut ), "=r" ( _addr )
                                        : "r" ( nb_loop ), "r" ( InterCoef ), "S" ( pSrc ), "0" ( pDestDebug ), "r" ( Delay32 ) //inputs
                                        : "%r15", "%r14", "%r8", "%r9", "%r10", "%xmm0", "%xmm1"
                                    );
                            RecoLutInd += nb_loop+1;
                        }// end ASM method
                        else //Intrinsics method*/
                        #endif //end #ifdef __ASMGCC64_OK
                        {
                            ////////////////////////////////////////////////////////////////////////////
                            //loop on the Slices for this region
                            for ( Slice = ( 2*SliceStart + j*2*NSlices ); Slice < ( 2*NSlices + j*2*NSlices ); Slice += 2 )
                            {
                                Time32 = RecoLut[RecoLutInd++];          // get time value in lambda*32
                                Time32 = Time32 + Delay32;               // Time32 corrected by delay
                                IndCoef = ( Time32 & 0x0006 ) + RecoLut[RecoLutInd++];
                                //Time32=(Time32>>3)&
                                data = ( __m128i * ) &ChannelStartPtr[Time32 >> 3];
                                //dataprefetch= (char *) (data + DeltaPrefetch);
                                //_mm_prefetch(dataprefetch,  _MM_HINT_NTA);
                                IQSum128i[Slice  ]=_mm_add_epi32 ( _mm_madd_epi16 ( _mm_loadu_si128 ( data ),InterCoef[IndCoef  ] ), IQSum128i[Slice  ] );
                                IQSum128i[Slice+1]=_mm_add_epi32 ( _mm_madd_epi16 ( _mm_loadu_si128 ( data ),InterCoef[IndCoef+1] ), IQSum128i[Slice+1] );
                            }   // End loop on slices
                        }
                    }       // End loop on Lines
                }           // End loop on Channels

                break;
        }               // End switch on modes
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        ////////////////////////////////////////////////////////////////////////////////////////////////////
// Compute and store output
        
        if( dataFormat == E_MATLAB )
        {
            switch ( OutputType )
            {
                case 0:     // (I,Q) output
                case 13:    // (I,Q) output rectangular
                case 10: // // (I,Q) output rectangular : we need to accumulate later in steerOutput
                    for ( j=0; j<NLines; j++ )
                    {
                        for ( i=0; i<2*NPixLine; i+=8 )
                        {
                            // Convert and add I values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                            // Convert and add Q values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                            _mm_empty();
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose Q values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                            // Intensity normalisation
                            xmmf128 = _mm_rcp_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                            xmme128 = _mm_mul_ps ( xmme128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                            xmma128 = _mm_mul_ps ( xmma128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                            
                            // Not needed for matlab format
                            //Interleave back I and Q
                            //xmmb128 = _mm_unpacklo_ps ( xmme128, xmma128 );    //Two first points
                            //xmmc128 = _mm_unpackhi_ps ( xmme128, xmma128 );    //Two next points
                            // Store points

//                             _mm_storeu_ps ( &IQPix[2*NPixLine*NLines* ( NRec-StartRec ) + 2*NPixLine*j+i], xmmb128 );
//                             _mm_storeu_ps ( &IQPix[2*NPixLine*NLines* ( NRec-StartRec ) + 2*NPixLine*j+i+4], xmmc128 );
                            
                            _mm_storeu_ps ( &IQbufferR[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2], xmme128 );
                            _mm_storeu_ps ( &IQbufferC[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2], xmma128 );
                        }
                    }
                    break;
                
                case 3:    // (I,Q) output rectangular accumulation // data are not steered later and we can accumulate now
                    for ( j=0; j<NLines; j++ )
                    {
                        for ( i=0; i<2*NPixLine; i+=8 )
                        {
                            // Convert and add I values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                            // Convert and add Q values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                            _mm_empty();
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose Q values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                            // Intensity normalisation
                            xmmf128 = _mm_rcp_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                            xmme128 = _mm_mul_ps ( xmme128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                            xmma128 = _mm_mul_ps ( xmma128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                            
                            // not needed for matlab format
                            //Interleave back I and Q
                            //xmmb128 = _mm_unpacklo_ps ( xmme128, xmma128 );    //Two first points
                            //xmmc128 = _mm_unpackhi_ps ( xmme128, xmma128 );    //Two next points
                            // Store points

//                             xmmg128 = _mm_loadu_ps(&IQPix[2*NPixLine*NLines* ( NRec-StartRec ) + 2*NPixLine*j+i]);
//                             xmmh128 = _mm_loadu_ps(&IQPix[2*NPixLine*NLines* ( NRec-StartRec ) + 2*NPixLine*j+i+4]);
                            
                            xmmg128 = _mm_loadu_ps(&IQPix[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2]);
                            xmmh128 = _mm_loadu_ps(&IQPixImag[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2]);

//                             xmmg128 = _mm_add_ps(xmmg128,xmmb128);
//                             xmmh128 = _mm_add_ps(xmmh128,xmmc128);
                            
                            xmmg128 = _mm_add_ps(xmmg128,xmme128);
                            xmmh128 = _mm_add_ps(xmmh128,xmma128);


                            _mm_storeu_ps ( &IQPix[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2], xmmg128 );
                            _mm_storeu_ps ( &IQPixImag[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2], xmmh128 );
                        }
                    }
                    break;
                case 1 :  //Intensity type
                    for ( j=0; j<NLines; j++ )
                    {
                        for ( i=0; i<2*NPixLine; i+=8 )
                        {
                            // Convert and add I values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                            // Convert and add Q values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                            _mm_empty();
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                            // Intensity computation
                            xmmb128 = _mm_sqrt_ps ( _mm_add_ps ( _mm_mul_ps ( xmma128, xmma128 ),_mm_mul_ps ( xmme128, xmme128 ) ) );
                            // Intensity normalisation
                            xmmf128 = _mm_rcp_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                            xmmb128 = _mm_mul_ps ( xmmb128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                            // Store result
                            _mm_storeu_ps ( &IQPix[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2], xmmb128 );
                        }
                    }
                    break;
                case 2 :    //Output square root of intensity on one short
                    for ( j=0; j<NLines; j++ )
                    {
                        for ( i=0; i<2*NPixLine; i+=8 )
                        {
                            // Convert and add I values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                            // Convert and add Q values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                            _mm_empty();
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                            // Intensity computation
                            xmmc128 = _mm_sqrt_ps ( _mm_add_ps ( _mm_mul_ps ( xmma128, xmma128 ),_mm_mul_ps ( xmme128, xmme128 ) ) );
                            // Intensity normalisation
                            xmmf128 = _mm_rcp_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                            xmmc128 = _mm_mul_ps ( xmmc128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                            // Square root of intensity computation
                            xmmb128 = _mm_min_ps ( _mm_sqrt_ps ( xmmc128 ), _mm_set1_ps ( 65535 ) );
                            //xmmb128 = _mm_sqrt_ps(xmmc128);
                            xmma128i = _mm_cvtps_epi32 ( xmmb128 );
                            IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2  ] = ( short ) _mm_extract_epi16 ( xmma128i, 0 );
                            IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2+1] = ( short ) _mm_extract_epi16 ( xmma128i, 2 );
                            IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2+2] = ( short ) _mm_extract_epi16 ( xmma128i, 4 );
                            IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2+3] = ( short ) _mm_extract_epi16 ( xmma128i, 6 );
                        }
                    }
                    _mm_empty();
                    break;
                    ///////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////////
                    // Output types with square root normalization
                    /////////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////////

                case 100:     // (I,Q) output
                case 110:    // (I,Q) output rectangular
                    for ( j=0; j<NLines; j++ )
                    {
                        for ( i=0; i<2*NPixLine; i+=8 )
                        {
                            // Convert and add I values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                            // Convert and add Q values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                            _mm_empty();
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose Q values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                            // Intensity normalisation
    //                             xmmf128 = _mm_rcp_ps(_mm_add_ps(_mm_min_ps(AntennaGain128[3*i/8], MaxRightGain128), _mm_min_ps(AntennaGain128[3*i/8+1], MaxLeftGain128)));
                            xmmf128 = _mm_rsqrt_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                            xmme128 = _mm_mul_ps ( xmme128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                            xmma128 = _mm_mul_ps ( xmma128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                            
                            // not needed in matlab format
                            //Interleave back I and Q
                            //xmmb128 = _mm_unpacklo_ps ( xmme128, xmma128 );    //Two first points
                            //xmmc128 = _mm_unpackhi_ps ( xmme128, xmma128 );    //Two next points
                            // Store points
//                             _mm_storeu_ps ( &IQPix[2*NPixLine*NLines* ( NRec-StartRec ) + 2*NPixLine*j+i], xmmb128 );
//                             _mm_storeu_ps ( &IQPix[2*NPixLine*NLines* ( NRec-StartRec ) + 2*NPixLine*j+i+4], xmmc128 );
                            
                            _mm_storeu_ps ( &IQPix[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2], xmme128 );
                            _mm_storeu_ps ( &IQPixImag[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2], xmma128 );
                        }
                    }
                    break;
                case 101 :  //Intensity type
                    for ( j=0; j<NLines; j++ )
                    {
                        for ( i=0; i<2*NPixLine; i+=8 )
                        {
                            // Convert and add I values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                            // Convert and add Q values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                            _mm_empty();
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                            // Intensity computation
                            xmmb128 = _mm_sqrt_ps ( _mm_add_ps ( _mm_mul_ps ( xmma128, xmma128 ),_mm_mul_ps ( xmme128, xmme128 ) ) );
                            // Intensity normalisation
                            xmmf128 = _mm_rsqrt_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                            xmmb128 = _mm_mul_ps ( xmmb128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                            // Store result
                            _mm_storeu_ps ( &IQPix[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2], xmmb128 );
                        }
                    }
                    break;
                case 102 :    //Output square root of intensity on one short
                    for ( j=0; j<NLines; j++ )
                    {
                        for ( i=0; i<2*NPixLine; i+=8 )
                        {
                            // Convert and add I values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                            // Convert and add Q values for 4 pixels
                            xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                            xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                            xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                            xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                            _mm_empty();
                            _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                            xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                            xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                            // Intensity computation
                            xmmc128 = _mm_sqrt_ps ( _mm_add_ps ( _mm_mul_ps ( xmma128, xmma128 ),_mm_mul_ps ( xmme128, xmme128 ) ) );
                            // Intensity normalisation
                            xmmf128 = _mm_rsqrt_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                            xmmc128 = _mm_mul_ps ( xmmc128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                            // Square root of intensity computation
                            xmmb128 = _mm_min_ps ( _mm_sqrt_ps ( xmmc128 ), _mm_set1_ps ( 65535 ) );
                            //xmmb128 = _mm_sqrt_ps(xmmc128);
                            xmma128i = _mm_cvtps_epi32 ( xmmb128 );
                            IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2  ] = ( short ) _mm_extract_epi16 ( xmma128i, 0 );
                            IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2+1] = ( short ) _mm_extract_epi16 ( xmma128i, 2 );
                            IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2+2] = ( short ) _mm_extract_epi16 ( xmma128i, 4 );
                            IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2+3] = ( short ) _mm_extract_epi16 ( xmma128i, 6 );
                        }
                    }
                    _mm_empty();
                    break;
            }
        
        }
        else
        {
            switch ( OutputType )
            {
            case 0:     // (I,Q) output
            case 10:    // (I,Q) output rectangular
                for ( j=0; j<NLines; j++ )
                {
                    for ( i=0; i<2*NPixLine; i+=8 )
                    {
                        // Convert and add I values for 4 pixels
                        xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                        xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                        xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                        xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                        _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                        xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                        xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                        // Convert and add Q values for 4 pixels
                        xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                        xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                        xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                        xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                        _mm_empty();
                        _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose Q values before addition
                        xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                        // Intensity normalisation
                        xmmf128 = _mm_rcp_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                        xmme128 = _mm_mul_ps ( xmme128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                        xmma128 = _mm_mul_ps ( xmma128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                        //Interleave back I and Q
                        xmmb128 = _mm_unpacklo_ps ( xmme128, xmma128 );    //Two first points
                        xmmc128 = _mm_unpackhi_ps ( xmme128, xmma128 );    //Two next points
                        // Store points
                        
                        _mm_storeu_ps ( &IQPix[2*NPixLine*NLines* ( NRec-StartRec ) + 2*NPixLine*j+i], xmmb128 );
                        _mm_storeu_ps ( &IQPix[2*NPixLine*NLines* ( NRec-StartRec ) + 2*NPixLine*j+i+4], xmmc128 );
                    }
                }
                break;
            case 1 :  //Intensity type
                for ( j=0; j<NLines; j++ )
                {
                    for ( i=0; i<2*NPixLine; i+=8 )
                    {
                        // Convert and add I values for 4 pixels
                        xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                        xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                        xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                        xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                        _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                        xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                        xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                        // Convert and add Q values for 4 pixels
                        xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                        xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                        xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                        xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                        _mm_empty();
                        _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                        xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                        // Intensity computation
                        xmmb128 = _mm_sqrt_ps ( _mm_add_ps ( _mm_mul_ps ( xmma128, xmma128 ),_mm_mul_ps ( xmme128, xmme128 ) ) );
                        // Intensity normalisation
                        xmmf128 = _mm_rcp_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                        xmmb128 = _mm_mul_ps ( xmmb128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                        // Store result
                        _mm_storeu_ps ( &IQPix[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2], xmmb128 );
                    }
                }
                break;
            case 2 :    //Output square root of intensity on one short
                for ( j=0; j<NLines; j++ )
                {
                    for ( i=0; i<2*NPixLine; i+=8 )
                    {
                        // Convert and add I values for 4 pixels
                        xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                        xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                        xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                        xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                        _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                        xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                        xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                        // Convert and add Q values for 4 pixels
                        xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                        xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                        xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                        xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                        _mm_empty();
                        _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                        xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                        // Intensity computation
                        xmmc128 = _mm_sqrt_ps ( _mm_add_ps ( _mm_mul_ps ( xmma128, xmma128 ),_mm_mul_ps ( xmme128, xmme128 ) ) );
                        // Intensity normalisation
                        xmmf128 = _mm_rcp_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                        xmmc128 = _mm_mul_ps ( xmmc128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                        // Square root of intensity computation
                        xmmb128 = _mm_min_ps ( _mm_sqrt_ps ( xmmc128 ), _mm_set1_ps ( 65535 ) );
                        //xmmb128 = _mm_sqrt_ps(xmmc128);
                        xmma128i = _mm_cvtps_epi32 ( xmmb128 );
                        IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2  ] = ( short ) _mm_extract_epi16 ( xmma128i, 0 );
                        IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2+1] = ( short ) _mm_extract_epi16 ( xmma128i, 2 );
                        IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2+2] = ( short ) _mm_extract_epi16 ( xmma128i, 4 );
                        IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2+3] = ( short ) _mm_extract_epi16 ( xmma128i, 6 );
                    }
                }
                _mm_empty();
                break;
                ///////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////
                // Output types with square root normalization
                /////////////////////////////////////////////////////////////////////////////////////
                ////////////////////////////////////////////////////////////////////////////////////

            case 100:     // (I,Q) output
            case 110:    // (I,Q) output rectangular
                for ( j=0; j<NLines; j++ )
                {
                    for ( i=0; i<2*NPixLine; i+=8 )
                    {
                        // Convert and add I values for 4 pixels
                        xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                        xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                        xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                        xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                        _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                        xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                        xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                        // Convert and add Q values for 4 pixels
                        xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                        xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                        xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                        xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                        _mm_empty();
                        _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose Q values before addition
                        xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                        // Intensity normalisation
//                             xmmf128 = _mm_rcp_ps(_mm_add_ps(_mm_min_ps(AntennaGain128[3*i/8], MaxRightGain128), _mm_min_ps(AntennaGain128[3*i/8+1], MaxLeftGain128)));
                        xmmf128 = _mm_rsqrt_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                        xmme128 = _mm_mul_ps ( xmme128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                        xmma128 = _mm_mul_ps ( xmma128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                        //Interleave back I and Q
                        xmmb128 = _mm_unpacklo_ps ( xmme128, xmma128 );    //Two first points
                        xmmc128 = _mm_unpackhi_ps ( xmme128, xmma128 );    //Two next points
                        // Store points
                        _mm_storeu_ps ( &IQPix[2*NPixLine*NLines* ( NRec-StartRec ) + 2*NPixLine*j+i], xmmb128 );
                        _mm_storeu_ps ( &IQPix[2*NPixLine*NLines* ( NRec-StartRec ) + 2*NPixLine*j+i+4], xmmc128 );
                    }
                }
                break;
            case 101 :  //Intensity type
                for ( j=0; j<NLines; j++ )
                {
                    for ( i=0; i<2*NPixLine; i+=8 )
                    {
                        // Convert and add I values for 4 pixels
                        xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                        xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                        xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                        xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                        _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                        xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                        xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                        // Convert and add Q values for 4 pixels
                        xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                        xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                        xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                        xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                        _mm_empty();
                        _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                        xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                        // Intensity computation
                        xmmb128 = _mm_sqrt_ps ( _mm_add_ps ( _mm_mul_ps ( xmma128, xmma128 ),_mm_mul_ps ( xmme128, xmme128 ) ) );
                        // Intensity normalisation
                        xmmf128 = _mm_rsqrt_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                        xmmb128 = _mm_mul_ps ( xmmb128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                        // Store result
                        _mm_storeu_ps ( &IQPix[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2], xmmb128 );
                    }
                }
                break;
            case 102 :    //Output square root of intensity on one short
                for ( j=0; j<NLines; j++ )
                {
                    for ( i=0; i<2*NPixLine; i+=8 )
                    {
                        // Convert and add I values for 4 pixels
                        xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 0], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 1] ); // Convert I values for 1st pixel to float
                        xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 4], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 5] ); // Convert I values for 1st pixel to float
                        xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 8], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 9] ); // Convert I values for 1st pixel to float
                        xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +12], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +13] ); // Convert I values for 1st pixel to float
                        _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                        xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                        xmme128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for I
                        // Convert and add Q values for 4 pixels
                        xmma128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 2], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 3] ); // Convert Q values for 1st pixel to float
                        xmmb128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 6], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) + 7] ); // Convert Q values for 1st pixel to float
                        xmmc128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +10], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +11] ); // Convert Q values for 1st pixel to float
                        xmmd128 = _mm_cvtpi32x2_ps ( IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +14], IQSum64[2* ( j*2*NSlices*NPixSlice + i ) +15] ); // Convert Q values for 1st pixel to float
                        _mm_empty();
                        _MM_TRANSPOSE4_PS ( xmma128, xmmb128, xmmc128, xmmd128 );        //transpose I values before addition
                        xmma128 = _mm_add_ps ( xmma128, xmmb128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmc128 );
                        xmma128 = _mm_add_ps ( xmma128, xmmd128 );                               // Final addition for Q
                        // Intensity computation
                        xmmc128 = _mm_sqrt_ps ( _mm_add_ps ( _mm_mul_ps ( xmma128, xmma128 ),_mm_mul_ps ( xmme128, xmme128 ) ) );
                        // Intensity normalisation
                        xmmf128 = _mm_rsqrt_ps ( _mm_add_ps ( _mm_min_ps ( AntennaGain128[3*i/8], MaxRightGain128 ), _mm_min_ps ( AntennaGain128[3*i/8+1], MaxLeftGain128 ) ) );
                        xmmc128 = _mm_mul_ps ( xmmc128, _mm_mul_ps ( xmmf128, AntennaGain128[3*i/8+2] ) );
                        // Square root of intensity computation
                        xmmb128 = _mm_min_ps ( _mm_sqrt_ps ( xmmc128 ), _mm_set1_ps ( 65535 ) );
                        //xmmb128 = _mm_sqrt_ps(xmmc128);
                        xmma128i = _mm_cvtps_epi32 ( xmmb128 );
                        IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2  ] = ( short ) _mm_extract_epi16 ( xmma128i, 0 );
                        IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2+1] = ( short ) _mm_extract_epi16 ( xmma128i, 2 );
                        IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2+2] = ( short ) _mm_extract_epi16 ( xmma128i, 4 );
                        IQPixshort[NPixLine*NLines* ( NRec-StartRec ) + NPixLine*j+i/2+3] = ( short ) _mm_extract_epi16 ( xmma128i, 6 );
                    }
                }
                _mm_empty();
                break;
            }
            
        }
        
    }
    //         Steer the data in array if rect output type.
    
    if( dataFormat == E_MATLAB )
    {
    
        if ( 10 == ( OutputType%100 ) || OutputType ==13 )
        {
            Lambda=SoundSpeed/DemodFrequ/1000;                      // Compute Lambda in mm
            JZero = int ( 65536*ZOrigin*tan ( SteerAngle ) / ( PiezoPitch*LinePitch ) +0.5 );   //Column abscissa at Z origine
            DeltaJ = int ( 65536*PitchPix*Lambda*tan ( SteerAngle ) / ( PiezoPitch*LinePitch ) +0.5 );     //DeltaJ for DeltaZ = 1
            if ( 0 != DeltaJ )
            {
              //  mexPrintf("SteerAngle = %g\n",SteerAngle);
                // accum if OutputType == 13
                SteerOutputFloat (IQbufferR, IQPix, NPixLine, NRecon*NLines, JZero, DeltaJ, ( int * ) &IQSum128i[0] , OutputType ==13);  // Steer data in the output array
                SteerOutputFloat (IQbufferC, IQPixImag, NPixLine, NRecon*NLines, JZero, DeltaJ, ( int * ) &IQSum128i[0], OutputType ==13);  // Steer data in the output array
            }
        }
        return 0;
    }
    else
    {
        if ( 10 == ( OutputType%100 ) )
        {
            Lambda=SoundSpeed/DemodFrequ/1000;                      // Compute Lambda in mm
            JZero = int ( 65536*ZOrigin*tan ( SteerAngle ) / ( PiezoPitch*LinePitch ) +0.5 );   //Column abscissa at Z origine
            DeltaJ = int ( 65536*PitchPix*Lambda*tan ( SteerAngle ) / ( PiezoPitch*LinePitch ) +0.5 );     //DeltaJ for DeltaZ = 1
        
            if ( 0 != DeltaJ )
            {
                SteerOutput ( ( double * ) &IQPix[0], NPixLine, NRecon*NLines, JZero, DeltaJ, ( int * ) &IQSum128i[0] );  // Steer data in the output array
            }
        }
        return 0;
    }

}                           //end of function
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
// function to steer output data in a rectangular array
// this is used for synthetic aperture in flat mode
void sp::SteerOutput ( double * IQArray, int NPixLine, int NCols, int JZero, int DeltaJ, int * Offset )
{
    int i, j, jrect, ind;
    int * IBottom;
    IBottom = ( int * ) &Offset[NPixLine];          // Pointer to last pixel on column.
    // Initialization IBottom array used for deepest pixel on every column
    jrect = ( JZero+32768 ) >>16;
    if ( DeltaJ > 0 )                                // case of positive steering angle
    {
        if ( jrect > NCols ) jrect = NCols;
        for ( j=0; j<jrect; j++ )
        {
            IBottom[j]= 0 ;
        }
        for ( j=jrect; j<NCols; j++ )
        {
            IBottom[j]= NPixLine ;
        }
    }
    if ( DeltaJ <  0 )                                // case of negative steering angle
    {
        jrect=NCols-1+jrect;
        if ( jrect < 0 ) jrect = -1;
        for ( j=0; j<=jrect; j++ )
        {
            IBottom[j]= NPixLine ;
        }
        for ( j=jrect+1; j<NCols; j++ )
        {
            IBottom[j]= 0 ;
        }
    }
    // Then we compute offsets for each line
    for ( i=0; i<NPixLine; i++ )
    {
        jrect = ( ( JZero+DeltaJ*i+32768 ) >>16 );
        if ( ( jrect<NCols ) && ( jrect> -NCols ) )
        {
            Offset[i] = -jrect*NPixLine;
            if ( ( jrect <= 0 ) && ( DeltaJ<0 ) ) jrect=NCols-1+jrect;
            IBottom[jrect] = i+1;
        }
    }
    // then we move data
    if ( DeltaJ <  0 )
    {
        for ( j=0; j<NCols; j++ )
        {
            ind = j*NPixLine;
            for ( i=0; i<IBottom[j]; i++ )
            {
                IQArray[ind] = IQArray[Offset[i] + ind];
                ind++;
            }
            for ( i=IBottom[j]; i<NPixLine; i++ )
            {
                IQArray[ind++] = 0;
            }
        }
    }
    if ( DeltaJ >  0 )
    {
        for ( j=NCols-1; j>=0; j-- )
        {
            ind = j*NPixLine;
            for ( i=0; i<IBottom[j]; i++ )
            {
                IQArray[ind] = IQArray[Offset[i] + ind];
                ind++;
            }
            for ( i=IBottom[j]; i<NPixLine; i++ )
            {
                IQArray[ind++] = 0;
            }
        }
    }
}       // end function SteerOutput

void sp::SteerOutputFloat ( float * IQArrayIn, float * IQArrayOut, int NPixLine, int NCols, int JZero, int DeltaJ, int * Offset , bool accum)
{
    int i, j, jrect, ind;
    int * IBottom;
    IBottom = ( int * ) &Offset[NPixLine];          // Pointer to last pixel on column.
    // Initialization IBottom array used for deepest pixel on every column
    jrect = ( JZero+32768 ) >>16;
    if ( DeltaJ > 0 )                                // case of positive steering angle
    {
        if ( jrect > NCols ) jrect = NCols;
        for ( j=0; j<jrect; j++ )
        {
            IBottom[j]= 0 ;
        }
        for ( j=jrect; j<NCols; j++ )
        {
            IBottom[j]= NPixLine ;
        }
    }
    if ( DeltaJ <  0 )                                // case of negative steering angle
    {
        jrect=NCols-1+jrect;
        if ( jrect < 0 ) jrect = -1;
        for ( j=0; j<=jrect; j++ )
        {
            IBottom[j]= NPixLine ;
        }
        for ( j=jrect+1; j<NCols; j++ )
        {
            IBottom[j]= 0 ;
        }
    }
    // Then we compute offsets for each line
    for ( i=0; i<NPixLine; i++ )
    {
        jrect = ( ( JZero+DeltaJ*i+32768 ) >>16 );
        if ( ( jrect<NCols ) && ( jrect> -NCols ) )
        {
            Offset[i] = -jrect*NPixLine;
            if ( ( jrect <= 0 ) && ( DeltaJ<0 ) ) jrect=NCols-1+jrect;
            IBottom[jrect] = i+1;
        }
    }
    // then we move data or accum
    
    if(accum)
    {
        if ( DeltaJ <  0 )
        {
            for ( j=0; j<NCols; j++ )
            {
                ind = j*NPixLine;
                for ( i=0; i<IBottom[j]; i++ )
                {
                    IQArrayOut[ind] += IQArrayIn[Offset[i] + ind++];
                }
            }
        }
        if ( DeltaJ >  0 )
        {
            for ( j=NCols-1; j>=0; j-- )
            {
                ind = j*NPixLine;
                for ( i=0; i<IBottom[j]; i++ )
                {
                    IQArrayOut[ind] += IQArrayIn[Offset[i] + ind++];
                }
            }
        }

    }
    else
    {
        if ( DeltaJ <  0 )
        {
            for ( j=0; j<NCols; j++ )
            {
                ind = j*NPixLine;
                for ( i=0; i<IBottom[j]; i++ )
                {
                    IQArrayOut[ind] = IQArrayIn[Offset[i] + ind++];
                }
                for ( i=IBottom[j]; i<NPixLine; i++ )
                {
                    IQArrayOut[ind++] = 0;
                }
            }
        }
        if ( DeltaJ >  0 )
        {
            for ( j=NCols-1; j>=0; j-- )
            {
                ind = j*NPixLine;
                for ( i=0; i<IBottom[j]; i++ )
                {
                    IQArrayOut[ind] = IQArrayIn[Offset[i] + ind++];
                }
                for ( i=IBottom[j]; i<NPixLine; i++ )
                {
                    IQArrayOut[ind++] = 0;
                }
            }
        }
    }
}       // end function SteerOutput
