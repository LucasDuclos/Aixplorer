
#include "CompLutSplitted.h"

int CompLutLinearPhotoacoustics(double *Probe, double *AcqInfo, double *ReconInfo, int *ReconMux, int *ReconAbsc, double * ReconDelay, short *Lut)
{
    int NPiezos, NLines, NPixLine, NSlices, NPixSlice, NRecon, PlaneWave, BandWidth, TxAlongRx, Folding, ImmediateLut;
    int FirstSample, NSamples, NChannels, MaxChannelsOneSide, PipeOffset, LutSizeBytes, NeededLutSizeBytes;
    double PiezoWidth, PiezoPitch, XOrigin, LinePitch, SteerAngle,RcvAngle, DemodFrequ, SoundSpeed;
    double ZOrigin, PitchPix, BandInterp, FootInterp, PeakDelay, Delay, MaxDelay, MinDelay,RXApodOne, RXApodZero;
    int Phase, PhaseI, PhaseQ , Time, Apod, Ech;
    const double Pi=M_PI;
    double x, Nphases, Filtre4[64*4], SumFiltre, Sin, Cos;
    double RXApod[4096], AvSoundSpeed[1024], *ProbeOut, *AcqInfoOut, *ReconInfoOut;
    float * AntennaGain;
    int *Lutint, *GlobalLut, *RecoLutMoinsStart, GlobalLutSize, InterCoefSize, IQSumSize, RecoLutSize;
    int ProbeSize, AcqInfoSize, ReconInfoSize, ReconInfoNext, OutputType, AntennaGainSize;
    short *InterCoef, *RecoLut;
    int i, NRec, Mode, ModeOk, Ardi1, Ardi2, CentralChannel, SyntAcq, ChannelOffset, Accum, ForcedMode;
    int AvailableModes[] = {0x1110308, 0x1210308, 0x1010308, 0x0010308, 0x0210308, 0x1210504, 0x1010504, 0x0010504, 0x0210504, 0x0010100,
                            0x1120504, 0x1220504, 0x1020504, 0x0020504, 0x0220504, 0x0020100, 0};         // Definition of allowed modes

    // this is the accumulation strategy for accumulated planewave
    int PlaneWaveModes[]    =  {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 ,-1};     // allowed accumulation choices
    int AccumPackets[]      =  {0, 0, 1, 1, 1, 1, 1, 1, 1, 1,  2,  1,  2,  1,  2,  3,  2 };     // number of accumulation packets
    int PostDivideBy2[]     =  {0, 0, 0, 0, 0, 1, 1, 1, 1, 1,  0,  1,  0,  1,  0,  0,  1 };     // final shift right (divide by 2)
    float InterpolationGain;
    __m128i *InterCoef128i, *RecoLut128i;
    __m128 xmma128, xmmb128, *AntennaGain128;

    #define MAX_INTERPOLATION_GAIN 256             //Interpolation coef limited to +/- 10 bits
    #define MaxGain 32.0                       //Gain multiplier for output
//     #define FFactorFlat 1.0                     // Coefficient for weighting in flat

    ////////////////////////////////////////////////////////////////////////////////////////////////
    /* Get Structures parameters */
    Lutint = (int *) &Lut[0];   // integer pointer to start of Lut
    // 'Probe'
    NPiezos = (int)Probe[0];    // integer unit. number of piezos on probe
    PiezoWidth = Probe[1];      // mm unit. Piezo width.
    PiezoPitch = Probe[2];      // mm unit. Piezo pitch.

    // 'AcqInfo'
    SteerAngle = AcqInfo[0] ;    // radian unit. Steering angle.
    BandWidth = (int)AcqInfo[1] ;     // 200 100 50. (I,Q) acquisition bandwidth.
        // FrequDivider not needed here
    DemodFrequ = AcqInfo[3] ;    // MHz unit. Demodulation frequency.
    SyntAcq = (int)AcqInfo[4];       // Synthetic acquisition flag. 0 or 1.
    ZOrigin = AcqInfo[5] ;       // mm unit. depth of first reconstructed points.
        // Depth not needed here
        // InitialTime not needed here
    FirstSample = (int)AcqInfo[8] ;   // integer. Index of first sample sent to PC (useful for partial transmission).
    NSamples = (int)AcqInfo[9] ;      // integer. Number of samples per channel.
        // NAcqChannels not needed here

    // 'ReconInfo
    NRecon = (int)ReconInfo[0] ;        // integer unit. Number of reconstructions
    PlaneWave = (int)(ReconInfo[1]-10);     // 0 or 1 . 1 if plane wave reconstruction, 0 if focalized reconstruction, 11 or 10 if photoacoustics
    NLines = (int)ReconInfo[2] ;   // integer unit. Number of lines in AcqInfo. 1 for plane wave, >= 1 for focused wave.
    XOrigin = ReconInfo[3] ;       // PiezoPitch unit. Left line abscissa offset/left edge of reference piezo. usually 0 for plane Wave.
    LinePitch = ReconInfo[4] ;     // PiezoPitch unit. Horizontal distance between two lines. 1 for normal resolution, 0.5 for double resolution.
    PitchPix = ReconInfo[5] ;      // Lambda unit. Distance between two pixels on reconstructed line.
    NPixLine = (int)ReconInfo[6] ;  // integer unit. Number of reconstructed slices
    PeakDelay = ReconInfo[7] ;          // µsec unit . Delay between start of transmit pulse, and maximum of transmit power.
    PipeOffset = (int)ReconInfo[8] ;    // integer. demod cycle unit . Time corresponding to first sample, assumed to be a 'I'.
    NChannels = (int)ReconInfo[9] ;    // integer. Number of acquisition channels.
    MaxChannelsOneSide = (int)ReconInfo[10] ;   // Maximum number of channels contributing on one side of a region.
    ChannelOffset = (int)ReconInfo[11]; // Channel offset in receive buffer.
    TxAlongRx = (int)ReconInfo[12];     // 1 = Receive aperture is set around firing abscissa.
    if (TxAlongRx >=10) {
        RcvAngle=0;
        TxAlongRx = TxAlongRx-10;
    }
    else {
        RcvAngle=SteerAngle;
    }
    BandInterp = ReconInfo[13] ;         // real [0 100]. Low pass filter for interpolation. 97 or 65
    FootInterp = ReconInfo[14] ;         // real . 100 * Foot of Hanning interpolation. 50 or 14
    RXApodOne = ReconInfo[15];          // max aperture value for full gain (Arditty coefficient)
    RXApodZero = ReconInfo[16];         // max aperture value for a contributing channel
    OutputType = (int)ReconInfo[17] ;   // Type of output : 0=(I,Q)   1= intensity
    NPixSlice = (int)ReconInfo[18];     // integer unit. Number of reconstructed pixels per slice.
    LutSizeBytes = (int)ReconInfo[19] ; // Lut size in bytes
    Mode = (int)ReconInfo[20] ;          // Will be recomputed later
    SoundSpeed = ReconInfo[21] ;    // m/sec unit. Speed of sound.

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Parameter check and status table

    if ( PiezoWidth > PiezoPitch ) return 1;   // PiezoWidth must be smaller than piezo pitch
    //if ( NPixLine != 4*(NPixLine/4) ) return 2;   // NPixLine (total number of pixels recontructed per line) must be a multiple of 4
    if ( ZOrigin < 0 ) return 3;   // ZOrigin must be positive (no pixel inside probe)
    if ( (BandWidth != 200)&&(BandWidth != 100)&&(BandWidth != 50) ) return 5;   // Bandwidth values : 200, 100, or 50
    if ( (BandInterp < 0) || (BandInterp > 100) ) return 6;   // BandInterp (interpolation band) must be between 0 and 100
    if ( MaxChannelsOneSide > NChannels ) return 7;   // MaxChannelsOneSide must be <= NChannels
    if ( (BandWidth == 200)&&(NPixSlice > 1)&& (8*PitchPix != int(8*PitchPix)) ) return 11;   // Multipixel slicing can only be performed when PitchPix is a multiple of lamba/8 in 200% mode, lambda/2 in 100% mode, lambda in 50% mode
    if ( (BandWidth == 100)&&(NPixSlice > 1)&& (2*PitchPix != int(2*PitchPix)) ) return 11;   // Multipixel slicing can only be performed when PitchPix is a multiple of lamba/8 in 200% mode, lambda/2 in 100% mode, lambda in 50% mode
    if ( (BandWidth ==  50)&&(NPixSlice > 1)&& (  PitchPix != int(  PitchPix)) ) return 11;   // Multipixel slicing can only be performed when PitchPix is a multiple of lamba/8 in 200% mode, lambda/2 in 100% mode, lambda in 50% mode
    if ( (OutputType != 0) && (OutputType != 1) && (OutputType != 2)&& (OutputType != 10)           // Output type must be 0,1, 2, 10, 100,101,102, or 110
            &&(OutputType != 100) && (OutputType != 101) && (OutputType != 102)&& (OutputType != 110)
            &&(OutputType != 200) && (OutputType != 201) && (OutputType != 202)&& (OutputType != 210)) return 4414;
    if ( ((OutputType%100) >= 10) && ((OutputType%100) < 20) && (TxAlongRx >= 10) ) return 23;

    // check for valid PlaneWave parameteravailable accumulation mode
    i=0;
    while ((PlaneWaveModes[i] > -1)&&(PlaneWaveModes[i] != PlaneWave))  {
        i=i+1;
    }
    if(PlaneWaveModes[i] == -1)  return 22;      // case where the parameter is invalid

    PlaneWave=PlaneWave + (AccumPackets[i]<<16) + (PostDivideBy2[i]<<24) ;  // pass parameters to reconstruct
//     printf("PlaneWave = %d\n", PlaneWave);
    // reduce interpolation gain if there is gain due to accumulation in planewavemode
    InterpolationGain =float(MAX_INTERPOLATION_GAIN);
//     printf("InterpolationGain0 = %f\n", InterpolationGain);
    //
    if (PlaneWave>1) {
        InterpolationGain=InterpolationGain/(PlaneWaveModes[i]/float(1<<PostDivideBy2[i]));
//         printf("PlaneWaveModes[i] = %d\n", PlaneWaveModes[i]);
//         printf("1<<PostDivideBy2[i] = %d\n", 1<<PostDivideBy2[i]);
        if (AccumPackets[i]>1)  InterpolationGain=InterpolationGain*2;
//         printf("InterpolationGain1 = %f\n", InterpolationGain);
    }


    //  0 : Lut Computation ok
    //  -1 : Mode not yet implemented
    //  1 : PiezoWidth must be smaller than PiezoPitch
    //  2 : available    NPixLine (total number of pixels recontructed per line) must be a multiple of 4
    //  3 : ZOrigin must be positive (no pixel inside probe)
    //  4 : RXApodOne and RXApodZero values must be between 0 and 2
    //  5 : Bandwidth values : 200, 100, or 50
    //  6 : BandInterp (interpolation band) must be between 0 and 100
    //  7 : MaxChannelsOneSide must be <= NChannels
    //  8 : ZOrigin too small for available RF data.
    //  9 : Image depth is too large for available RF data.
    // 10 : available
    // 11 : Multipixel slicing can only be performed when PitchPix is a multiple of lamba/8 in 200% mode, lambda/2 in 100% mode, lambda in 50% mode
    // 12 : Lut size is not large enough
    // 13 : RXApod must be >=0 and <= 1
    // 14 : Output type must be 0,1, 2, or 10
    // 15 : Region cannot be left of aperture
    // 16 : Region cannot be right of aperture
    // 17 : For rectangular output NPixSlice must be egal to 1
    // 18 : Inconsistent bandwidth for forced mode
    // 19 : Geometry does not allow folding
    // 20 : Variable delays do not allow an immediate LUT
    // 21 : NpixSlice not coherent with forced mode
    // 22 : Invalid PlaneWave value
    // 23 : Rectangular output type is not allowed when the output is already forced rectangular with TxAlongRx >=10
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Search for reconstruction mode Mode
    //////////////////////////
    // Mode format : ImmediateLut = 7th nib, Accum = bit 21, Folding= bit 20, BandWidth = 5th nib, NPixSlice = 3rd and 4th nib, 8*PithPix = 1st and 2nd nib
    /////////////////////////////////////////////////////////////////////////////////
    // First check if mode is imposed, and check also for special modes
    if (Mode > 0) {     //Case where mode is imposed
        ForcedMode = Mode;                               // Store desired mode
        Accum=(Mode & 0x00600000) >> 21;                  // Get Accum parameter on 2 bits
    }
    else if ( Mode==0) {
        ForcedMode = 0;  
        Accum=0;                                // An accumulation mode is not asked for
    }
    else if ( Mode<0) {
        ForcedMode = 0;  
        Accum=-Mode;                                // An accumulation mode is asked for
    }
    else {
        return -1;                              // Mode not yet implemented
    }
    //////////////////////////////////////////////////////////////////////////////////////
    // Now compute possible range of modes knowing reconstruction parameters
    // Set Bandwidth
    Mode = BandWidth/10;
    if (Mode >=10)
        Mode = Mode/10;                // Case of 200% and 100% modes
    Mode = (Mode << 16);                          // 5th nib is bandwidth
    // Set Accum bit
    Mode = Mode + (Accum << 21);
    // Check if geometry allows folding
    if ((NLines == 1) && ((XOrigin ==0)||(XOrigin == 0.5)) && (SteerAngle ==0)) {
        Folding =1;                             // It might be possible to perform folding
    }
    else {
        Folding =0;
    }
    // Check if firing delays allow the use of an immediate Lut
    // Compute maximum delay and check if delay is constant
    MaxDelay =0;
    MinDelay=100000;
    Delay = int(32*ReconDelay[0]*DemodFrequ+ 0.5);
    for (NRec = 0; NRec < NRecon; NRec++) {             //loop on all reconstructions
        if ( (ReconDelay[NRec]*DemodFrequ) > MaxDelay)  MaxDelay = ReconDelay[NRec]*DemodFrequ;
        if ( (ReconDelay[NRec]*DemodFrequ) < MinDelay)  MinDelay = ReconDelay[NRec]*DemodFrequ;
        if ( int(32*ReconDelay[NRec]*DemodFrequ + 0.5) != Delay)  Delay = - 100000;
    }
    if (Delay > -100000) {
        ImmediateLut = 1;                               // Delays are constant, an immediate Lut can be envisaged
    }
    else {
        ImmediateLut = 0;                                // Delays are not constant, an immediate Lut cannot be envisaged
    }
    // Check if rectangular output is asked (for synthetic mode)
    if ( ((OutputType%100) >= 10) && ((OutputType%100) < 20) ) {
       if (NPixSlice !=1 ) {
          return 17;
       }
    }
    ///////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////
    // Verify that in case of forced mode parameters are consistent with mode
    if (ForcedMode) {
        if ((ForcedMode & 0x000F0000)!=(Mode & 0x000F0000)) return 18; // error on bandwidth
        if (ForcedMode & 0x00100000) {          // case where we want to force folding
            if ((!Folding)) return 19;          // error because folding is not allowed
        }
        else Folding=0;                         // In that case we do not want folding
        if (ForcedMode & 0x0F000000) {         // case where we want to force an immediate Lut
            if ((!ImmediateLut)) return 20;     // error because Immediate Lut is not allowed
        }
        else ImmediateLut=0;                    // in that case we do not want an Immediate Lut
        if ( (NPixSlice==0) || (NPixSlice != (ForcedMode & 0x0000FF00)) >> 8) return 21;        // incoherent NPixslice
    }
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // Check mode availability for reconstruction
    i=0;
    ModeOk = 0;
    // Case where NPixSlice is imposed ; this may be a forced mode
    if (NPixSlice > 0) {
        NSlices = int(ceil(NPixLine/double(NPixSlice)));
        // Set NPixSlice
        if (NPixSlice > 255) return -1;              // Too much pixels per slice
        Mode = Mode + (NPixSlice << 8);             // 3rd and 4th nib for number of Pixels per slice
        // Set 8*Lambda for first and 2nd nib
        if (NPixSlice > 1)  {
            if (PitchPix >= 32)  return -1;          // Pixel pitch is too large
            Mode = Mode + int(8*PitchPix);          // 1st and 2nd digits are proportional to pixel pitch
        }
        while(AvailableModes[i] != 0) {
            // First search in case of forced mode
            if ((ForcedMode > 0)&&(ForcedMode == AvailableModes[i])&&(ModeOk ==0)) {
                // done in the case of forced mode
                ModeOk = 1;     // Recognized available mode
                Mode = AvailableModes[i];
            }
            // Now search in case of automatic mode
            else if (((Mode +(Folding<<20)+(ImmediateLut<<24)) == (AvailableModes[i]))&&(ModeOk ==0)) {
                // done if no immediate lut, or if small enough immediate lut
                ModeOk = 1;     // Recognized available mode
                Mode = AvailableModes[i];

            }
            else if (((Mode +(ImmediateLut<<24)) == (AvailableModes[i]))&&(ModeOk ==0)) {
                // done if no immediate lut, or if small enough immediate lut
                ModeOk = 1;     // Recognized available mode
                Mode = AvailableModes[i];
                Folding = 0;

            }
            else if (((Mode +(Folding<<20)) == (AvailableModes[i]))&&(ModeOk ==0)) {
                ModeOk = 1;     // Recognized available mode
                Mode = AvailableModes[i];
                ImmediateLut = 0;
            }
            else if ((Mode  == (AvailableModes[i]))&&(ModeOk ==0)) {
                ModeOk = 1;     // Recognized available mode
                Mode = AvailableModes[i];
                ImmediateLut = 0;
                Folding=0;
            }

            i++ ;
        }
    }
    // Case where NPixSlice is not imposed (NPixSlice = 0). Automatic NPixslice choice and automatic mode only
    else {
        // Case of possible multipixel slice.
        if ((8*PitchPix) == int(8*PitchPix)) {
            Mode = Mode + int(8*PitchPix);
            while(AvailableModes[i] != 0) {
                if (((Mode +(Folding<<20)+(ImmediateLut<<24)) == (AvailableModes[i]& 0x0FFF00FF))&&(ModeOk ==0)) {
                    NPixSlice = (AvailableModes[i] & 0x0000FF00) >> 8;
                    NSlices = int(ceil(NPixLine/double(NPixSlice)));
                    // done if no immediate lut, or if small enough immediate lut
                    ModeOk = 1;     // Recognized available mode
                    Mode = AvailableModes[i];
                }
                else if (((Mode +(ImmediateLut<<24)) == (AvailableModes[i]& 0x0FFF00FF))&&(ModeOk ==0)) {
                    NPixSlice = (AvailableModes[i] & 0x0000FF00) >> 8;
                    NSlices = int(ceil(NPixLine/double(NPixSlice)));
                    // done if no immediate lut, or if small enough immediate lut
                    ModeOk = 1;     // Recognized available mode
                    Mode = AvailableModes[i];
                    Folding = 0;
                }
                else if (((Mode +(Folding<<20)) == (AvailableModes[i]& 0x0FFF00FF))&&(ModeOk ==0)) {
                    NPixSlice = (AvailableModes[i] & 0x0000FF00) >> 8;
                    NSlices = int(ceil(NPixLine/double(NPixSlice)));
                    ModeOk = 1;     // Recognized available mode
                    Mode = AvailableModes[i];
                    ImmediateLut = 0;
                }
                else if ((Mode  == (AvailableModes[i]& 0x0FFF00FF))&&(ModeOk ==0)) {
                    NPixSlice = (AvailableModes[i] & 0x0000FF00) >> 8;
                    NSlices = int(ceil(NPixLine/double(NPixSlice)));
                    ModeOk = 1;     // Recognized available mode
                    Mode = AvailableModes[i];
                    ImmediateLut = 0;
                    Folding=0;
                }

                i++ ;
            }
        }
        // no multipixel slice possibility
        if (ModeOk ==0 ) {
            i=0;
            NPixSlice = 1;                          // Set NPixSlice and mode
            NSlices = int(ceil(NPixLine/double(NPixSlice)));
            Mode = (Mode & 0x002F0000) + (NPixSlice << 8);      //Set for NPixSlice = 1
            while(AvailableModes[i] != 0) {
                if (((Mode +(Folding<<20)+(ImmediateLut<<24)) == (AvailableModes[i]))&&(ModeOk ==0)) {
                    // done if no immediate lut, or if small enough immediate lut
                    ModeOk = 1;     // Recognized available mode
                    Mode = AvailableModes[i];
                }
                else if (((Mode +(ImmediateLut<<24)) == (AvailableModes[i]))&&(ModeOk ==0)) {
                    // done if no immediate lut, or if small enough immediate lut
                    ModeOk = 1;     // Recognized available mode
                    Mode = AvailableModes[i];
                    Folding = 0;
                }
                else if (((Mode +(Folding<<20)) == (AvailableModes[i]))&&(ModeOk ==0)) {
                    ModeOk = 1;     // Recognized available mode
                    Mode = AvailableModes[i];
                    ImmediateLut = 0;
                }
                else if ((Mode  == (AvailableModes[i]))&&(ModeOk ==0)) {
                    ModeOk = 1;     // Recognized available mode
                    Mode = AvailableModes[i];
                    ImmediateLut = 0;
                    Folding=0;
                }

                i++ ;
            }
        }
    }
    if (ModeOk != 1) return -1;

    // Set flag to be used later to compute correct antenna gain
    if ((Folding == 1) && (XOrigin == 0.5)) {
        CentralChannel=1;
    }
    else {
        CentralChannel=0;
    }
    // Allocate Structures space
    //Structures size in shorts
    ProbeSize = 3*4;
    AcqInfoSize = 11*4;
    ReconInfoSize = 26*4;
    ReconInfoNext=ReconInfoSize/4;
    if (BandWidth == 200) InterCoefSize = 2*4*16*8;
    if (BandWidth == 100) InterCoefSize = 2*16*8*8;
    IQSumSize=2*NSlices*NPixSlice*NLines*8;
    AntennaGainSize = 3*(NSlices*NPixSlice+4)*2;       // 1 float per pixel of center line on each side of line, plus one float for depth gain
    GlobalLutSize=NRecon*2;
//      RecoLutSize = (2-Folding)*(1+(2+ImmediateLut)*NSlices)*NLines*MaxChannelsOneSide*(1+7*ImmediateLut); // size /2 if folding, *8 if immediate lut
    RecoLutSize = (2-Folding)*(1+2*NSlices)*NLines*MaxChannelsOneSide; // size /2 if folding, *8 if immediate lut
    NeededLutSizeBytes = (ProbeSize + AcqInfoSize + ReconInfoSize + InterCoefSize + IQSumSize + GlobalLutSize + AntennaGainSize + 2 + RecoLutSize + 8)*2; // add 8 (one SSE) for RecoLut alignment
    if (LutSizeBytes < NeededLutSizeBytes) {
       Lutint[0] = NeededLutSizeBytes;
       return 12;                       // Lut size is not large enough
    }
    // Pointers to structures
    ProbeOut = (double *) &Lut[0];
    AcqInfoOut = (double *) &Lut[ProbeSize];
    ReconInfoOut = (double *) &Lut[ProbeSize + AcqInfoSize];
    InterCoef = (short *) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize];
    InterCoef128i = (__m128i *) &InterCoef[0];
//     GlobalLut = (int *) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize +InterCoefSize + IQSumSize];
//     AntennaGain = (float*) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize + InterCoefSize + IQSumSize + GlobalLutSize];
//     RecoLutMoinsStart = (int *) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize + InterCoefSize + IQSumSize + GlobalLutSize + AntennaGainSize];
    AntennaGain = (float*) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize +InterCoefSize + IQSumSize];
    AntennaGain128 =(__m128 *) &AntennaGain[0];
    GlobalLut = (int *) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize +InterCoefSize + IQSumSize + AntennaGainSize];
    RecoLutMoinsStart = (int *) &Lut[ProbeSize + AcqInfoSize + ReconInfoSize + InterCoefSize + IQSumSize + AntennaGainSize + GlobalLutSize];

    RecoLut = (short *) &Lut[((ProbeSize + AcqInfoSize + ReconInfoSize + InterCoefSize + IQSumSize + AntennaGainSize + GlobalLutSize + 2)& 0xFFFFFFF8) + 0x00000008];
    RecoLut128i = (__m128i *) &RecoLut[0];
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Compute RXApod array

    if(RXApodZero >0) {     // Case where trapezoidal shape is specified with two Arditty parameters. We must then create the apodization table
        Ardi1=int(1000*RXApodOne+0.5);
        Ardi2=int(1000*RXApodZero+0.5);
        if ( (Ardi1 < 0) || (Ardi1 > 2000) || (Ardi2 < 0) || (Ardi2 > 2000) ) return 4;   // RXApodOne and RXApodZero values must be between 0 and 2
        // Set maximum value
        for (i=0; i<Ardi2; i++) {
            RXApod[i] = 1;
        }
        // Case where there is effective apodization
        if (Ardi1<Ardi2) {
            for (i=Ardi1; i<Ardi2; i++) {
                RXApod[i]= 1- (i-Ardi1)/double(Ardi2-Ardi1);
            }
        }
        // Set gain to zero beyond max aperture
        for (i=Ardi2; i<2048; i++) {
            RXApod[i]= 0;
        }
        // Translate for positive values
        for (i=0; i<2048; i++) {
            RXApod[i+2048]= RXApod[i];
        }
        // Create negative values identical as positive
        RXApod[0]=0;
        for (i=1; i<2048; i++) {
            RXApod[i]= RXApod[4096-i];
        }
    }
    else
    {
        for (i=0; i<4096; i++) {
            RXApod[i]= ReconInfo[ReconInfoNext+i];
        }
        ReconInfoNext += 4096;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Compute InterCoef Array

    switch (BandWidth) {
        case 200:     // Compute InterCoef for 200 % bandwidth
            // First compute interpolation filter
            Nphases = 8;   // 8 interpolation states for 200% case
            for (Phase=0; Phase<Nphases; Phase++) {
                SumFiltre = 0;
                for (Ech=0; Ech<4; Ech++) {
                    x=Ech-1 - (0.5 + Phase)/Nphases ;
                    Filtre4[4*Phase+Ech]=sin(M_PI*x*BandInterp/100)/(M_PI*x*BandInterp/100)*(1+cos(M_PI/2*x)+FootInterp/50);
                    SumFiltre += Filtre4[4*Phase+Ech];
                }
                for (Ech=0; Ech<4; Ech++) {
                    Filtre4[4*Phase+Ech]= Filtre4[4*Phase+Ech]/SumFiltre;
                }
            }
            // Then combine with rotation, and with apodization
            for (Apod=0; Apod<16; Apod++) {         // compute only one fourth of coefficients
                for (Time=0; Time<4; Time++) {
                    PhaseI = Time;
                    PhaseQ = Time+4;
                    Sin = -sin(2*M_PI*Time/16);
                    Cos = cos(2*M_PI*Time/16);
                    // Interpolation for I
                    InterCoef[64*Apod+16*Time+0]= int(InterpolationGain*Filtre4[4*PhaseQ+0]*(-Sin)*( 1)*(Apod+0.5)/15.5 +0.5); //Q
                    InterCoef[64*Apod+16*Time+1]= int(InterpolationGain*Filtre4[4*PhaseI+0]*( Cos)*(-1)*(Apod+0.5)/15.5 +0.5); //-I
                    InterCoef[64*Apod+16*Time+2]= int(InterpolationGain*Filtre4[4*PhaseQ+1]*(-Sin)*(-1)*(Apod+0.5)/15.5 +0.5); //-Q
                    InterCoef[64*Apod+16*Time+3]= int(InterpolationGain*Filtre4[4*PhaseI+1]*( Cos)*( 1)*(Apod+0.5)/15.5 +0.5); //I
                    InterCoef[64*Apod+16*Time+4]= int(InterpolationGain*Filtre4[4*PhaseQ+2]*(-Sin)*( 1)*(Apod+0.5)/15.5 +0.5); //Q
                    InterCoef[64*Apod+16*Time+5]= int(InterpolationGain*Filtre4[4*PhaseI+2]*( Cos)*(-1)*(Apod+0.5)/15.5 +0.5); //-I
                    InterCoef[64*Apod+16*Time+6]= int(InterpolationGain*Filtre4[4*PhaseQ+3]*(-Sin)*(-1)*(Apod+0.5)/15.5 +0.5); //-Q
                    InterCoef[64*Apod+16*Time+7]= int(InterpolationGain*Filtre4[4*PhaseI+3]*( Cos)*( 1)*(Apod+0.5)/15.5 +0.5); //I
                    // Interpolation for Q
                    InterCoef[64*Apod+16*Time+8] = int(InterpolationGain*Filtre4[4*PhaseQ+0]*(Cos)*( 1)*(Apod+0.5)/15.5 +0.5); //Q
                    InterCoef[64*Apod+16*Time+9] = int(InterpolationGain*Filtre4[4*PhaseI+0]*(Sin)*(-1)*(Apod+0.5)/15.5 +0.5); //-I
                    InterCoef[64*Apod+16*Time+10]= int(InterpolationGain*Filtre4[4*PhaseQ+1]*(Cos)*(-1)*(Apod+0.5)/15.5 +0.5); //-Q
                    InterCoef[64*Apod+16*Time+11]= int(InterpolationGain*Filtre4[4*PhaseI+1]*(Sin)*( 1)*(Apod+0.5)/15.5 +0.5); //I
                    InterCoef[64*Apod+16*Time+12]= int(InterpolationGain*Filtre4[4*PhaseQ+2]*(Cos)*( 1)*(Apod+0.5)/15.5 +0.5); //Q
                    InterCoef[64*Apod+16*Time+13]= int(InterpolationGain*Filtre4[4*PhaseI+2]*(Sin)*(-1)*(Apod+0.5)/15.5 +0.5); //-I
                    InterCoef[64*Apod+16*Time+14]= int(InterpolationGain*Filtre4[4*PhaseQ+3]*(Cos)*(-1)*(Apod+0.5)/15.5 +0.5); //-Q
                    InterCoef[64*Apod+16*Time+15]= int(InterpolationGain*Filtre4[4*PhaseI+3]*(Sin)*( 1)*(Apod+0.5)/15.5 +0.5); //I
                }
            }
            break;
        ////////////////////////////////////////////////////////////////////////////////////////////
        case 100:       // Compute InterCoef for 100 % bandwidth
            // First compute interpolation filter
            Nphases = 16;   // 8 interpolation states for 200% case
            for (Phase=0; Phase<Nphases; Phase++) {
                SumFiltre = 0;
                for (Ech=0; Ech<4; Ech++) {
                    x=Ech-1 - (0.5 + Phase)/Nphases ;
                    Filtre4[4*Phase+Ech]=sin(M_PI*x*BandInterp/100)/(M_PI*x*BandInterp/100)*(1+cos(M_PI/2*x)+FootInterp/50);
                    SumFiltre += Filtre4[4*Phase+Ech];
                }
                for (Ech=0; Ech<4; Ech++) {
                    Filtre4[4*Phase+Ech]= Filtre4[4*Phase+Ech]/SumFiltre;
                }
            }
            // Then combine with rotation, and with apodization
            for (Apod=0; Apod<8; Apod++) {         // compute only one fourth of coefficients
                for (Time=0; Time<4; Time++) {
                    PhaseI = Time;
                    PhaseQ = Time+12;
                    Sin = -sin(2*M_PI*Time/16);
                    Cos = cos(2*M_PI*Time/16);
                    // Interpolation for I
                    InterCoef[256*Apod+16*Time+0]= int(InterpolationGain*Filtre4[4*PhaseQ+0]*(-Sin)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+1]= int(InterpolationGain*Filtre4[4*PhaseI+0]*( Cos)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+2]= int(InterpolationGain*Filtre4[4*PhaseQ+1]*(-Sin)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+3]= int(InterpolationGain*Filtre4[4*PhaseI+1]*( Cos)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+4]= int(InterpolationGain*Filtre4[4*PhaseQ+2]*(-Sin)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+5]= int(InterpolationGain*Filtre4[4*PhaseI+2]*( Cos)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+6]= int(InterpolationGain*Filtre4[4*PhaseQ+3]*(-Sin)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+7]= int(InterpolationGain*Filtre4[4*PhaseI+3]*( Cos)*(Apod+0.5)/7.5 +0.5); //I
                    // Interpolation for Q
                    InterCoef[256*Apod+16*Time+8] = int(InterpolationGain*Filtre4[4*PhaseQ+0]*(Cos)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+9] = int(InterpolationGain*Filtre4[4*PhaseI+0]*(Sin)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+10]= int(InterpolationGain*Filtre4[4*PhaseQ+1]*(Cos)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+11]= int(InterpolationGain*Filtre4[4*PhaseI+1]*(Sin)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+12]= int(InterpolationGain*Filtre4[4*PhaseQ+2]*(Cos)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+13]= int(InterpolationGain*Filtre4[4*PhaseI+2]*(Sin)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+14]= int(InterpolationGain*Filtre4[4*PhaseQ+3]*(Cos)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+15]= int(InterpolationGain*Filtre4[4*PhaseI+3]*(Sin)*(Apod+0.5)/7.5 +0.5); //I

                }
                for (Time=4; Time<16; Time++) {
                    PhaseI = Time;
                    PhaseQ = Time-4;
                    Sin = -sin(2*M_PI*Time/16);
                    Cos = cos(2*M_PI*Time/16);
                    // Interpolation for I
                    InterCoef[256*Apod+16*Time+0]= int(InterpolationGain*Filtre4[4*PhaseI+0]*( Cos)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+1]= int(InterpolationGain*Filtre4[4*PhaseQ+0]*(-Sin)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+2]= int(InterpolationGain*Filtre4[4*PhaseI+1]*( Cos)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+3]= int(InterpolationGain*Filtre4[4*PhaseQ+1]*(-Sin)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+4]= int(InterpolationGain*Filtre4[4*PhaseI+2]*( Cos)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+5]= int(InterpolationGain*Filtre4[4*PhaseQ+2]*(-Sin)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+6]= int(InterpolationGain*Filtre4[4*PhaseI+3]*( Cos)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+7]= int(InterpolationGain*Filtre4[4*PhaseQ+3]*(-Sin)*(Apod+0.5)/7.5 +0.5); //Q
                    // Interpolation for Q
                    InterCoef[256*Apod+16*Time+8] = int(InterpolationGain*Filtre4[4*PhaseI+0]*(Sin)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+9] = int(InterpolationGain*Filtre4[4*PhaseQ+0]*(Cos)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+10]= int(InterpolationGain*Filtre4[4*PhaseI+1]*(Sin)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+11]= int(InterpolationGain*Filtre4[4*PhaseQ+1]*(Cos)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+12]= int(InterpolationGain*Filtre4[4*PhaseI+2]*(Sin)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+13]= int(InterpolationGain*Filtre4[4*PhaseQ+2]*(Cos)*(Apod+0.5)/7.5 +0.5); //Q
                    InterCoef[256*Apod+16*Time+14]= int(InterpolationGain*Filtre4[4*PhaseI+3]*(Sin)*(Apod+0.5)/7.5 +0.5); //I
                    InterCoef[256*Apod+16*Time+15]= int(InterpolationGain*Filtre4[4*PhaseQ+3]*(Cos)*(Apod+0.5)/7.5 +0.5); //Q

                }
            }
            break;
        case 50:     // Compute InterCoef for 50 % bandwidth
            // to be done
            break;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Compute RecoLut : Indexed Lut

    int j, jrect, k, p, Indlut, Arditty, ApodInt, AntennaLeft, IndSample, IndCoef, DeltaJ=0, JZero=0;
    double SliceThickness, ZSliceCenter, XPiezo, XSlice, XLine, XfromBeamCenter;
//     float A, B;
    double Lambda, ReturnDistance, TransmitDistance;
    double TimeOfFlightCenterSlice, TimeOfFlightStartSlice, TimeCycle;
    short SliceStart;

    Indlut = 0;             // initialize LUT indice

    Lambda=SoundSpeed/DemodFrequ/1000;                      // Compute Lambda in mm
    if ( ((OutputType%100) >= 10) && ((OutputType%100) < 20) ) {        // Case of rectangular output
        JZero = int(65536*ZOrigin*tan(SteerAngle)/(PiezoPitch*LinePitch)+0.5);              //Column abscissa at Z origine
        DeltaJ = int(65536*PitchPix*Lambda*tan(SteerAngle)/(PiezoPitch*LinePitch)+0.5);     //DeltaJ for DeltaZ = 1
        SliceThickness = NPixSlice*PitchPix*Lambda;          // actually in that case NPixSlice = 1
    }
    else {
        SliceThickness = NPixSlice*PitchPix*Lambda*cos(RcvAngle);         // Depth of slice in mm.
    }

//     A = float(ZOrigin*MaxGain*MaxGain/(FFactorFlat*PiezoPitch*NChannels));                             // for computing aperture in receive
//     B = float(PitchPix*Lambda*cos(SteerAngle)*MaxGain*MaxGain/(FFactorFlat*PiezoPitch*NChannels));     // to compute aperture in receive
    //////////////////////////////////////////////////////////////
    // load real speed of sound if it exists. In that case SoundSpeed has a fractional part
    // This part is not yet completely debug, so just take constant sound speed for the moment
//     if(SoundSpeed > int(SoundSpeed)) {
//         for (i=0; i<1024; i++)  {
//             AvSoundSpeed[i]=ReconInfo[ReconInfoNext+i];
//         }
//         ReconInfoNext += 1024;
//     }
//     else {
//         for (i=0; i<1024; i++)  {
//             AvSoundSpeed[i]=SoundSpeed;
//         }
//     }
    for (i=0; i<1024; i++)  {
            AvSoundSpeed[i]=SoundSpeed;
        }
    ///////////////////////////////////////////////////////
    // Initialize AntennaGain array
    for (i=0; i<AntennaGainSize/2; i++) {               //loop on all pixels for central line on every side of line
        AntennaGain[i] = 0;
    }
    ////////////////////////////////////////////////////////
    // Start Recolut exploration
    for (k = 0; k<(2-Folding)*MaxChannelsOneSide; k ++) {         //Loop on NChannels on both sides of Region

        if ( k < MaxChannelsOneSide ) {                 // Case of piezo on right side of region mm units
            XPiezo = (k + 0.5)*PiezoPitch;              // Piezo abscissa mm units
            AntennaLeft = 0;                            // to accumulate Channel contribution on right
        }
        else {
            XPiezo = (MaxChannelsOneSide - k - 0.5)*PiezoPitch;     // Case of piezo on left side of region mm units
            AntennaLeft = 1;                            // to accumulate Channel contribution on the left
        }
        if ( k == MaxChannelsOneSide ) {                // Start of analysis on left side
            RecoLutMoinsStart[0] = Indlut;              // Offset for RecoLut on left side
        }

        for (j=0; j<NLines; j++) {                      // loop on all lines
            SliceStart = -1;                                 //Starting line slice to which contributes every channel
            for (i=0; i<NSlices; i++) {                      // loop on all slices
                ZSliceCenter = ZOrigin + SliceThickness*(i + ((double)(NPixSlice-1))/2/NPixSlice);
                XLine = (XOrigin + j*LinePitch)*PiezoPitch;     // Line origin Abscissa abscissa in mm
                if ( ((OutputType%100) >= 10) && ((OutputType%100) < 20) ) {        // Case of rectangular output
                    jrect = ((JZero+DeltaJ*i+32768)>>16) + j;                 // Column index closest to abscissa
                    XSlice = (XOrigin + jrect*LinePitch)*PiezoPitch;      // Abscissa on grid
                }
                else {
                    XSlice = XLine + ZSliceCenter*tan(RcvAngle);  // Slice center abscissa in mm
                }
                XfromBeamCenter = (2*j+1-NLines)*(LinePitch*PiezoPitch)/2;
                // Return distance in mm :
                ReturnDistance = sqrt((XPiezo-XSlice)*(XPiezo-XSlice) + ZSliceCenter*ZSliceCenter);
                Arditty=2048+int(1000*2*((XPiezo-XSlice + TxAlongRx*(XSlice - XLine))/(ReturnDistance+0.001))*PiezoWidth/Lambda + 0.5);    // Apodisation centered around line origin
                if (Arditty >= 4096) {
                    Arditty=4095;
                }
                if (Arditty < 0) {
                    Arditty=0;
                }
                if ((RXApod[Arditty] > 0) || (SliceStart>=0)) {          // Case where the channel contributes to the point
                    //TransmitDistance = ZSliceCenter/cos(SteerAngle)+XfromBeamCenter*sin(SteerAngle);
                    TransmitDistance = 0.0;// This is PAT and not pulse echo

//                     TimeOfFlightCenterSlice = (TransmitDistance + ReturnDistance)/Lambda;
                    TimeOfFlightCenterSlice = (TransmitDistance + ReturnDistance)*1000*DemodFrequ/AvSoundSpeed[int(TransmitDistance + ReturnDistance+0.5)];
                    TimeOfFlightStartSlice = TimeOfFlightCenterSlice - (NPixSlice-1)*PitchPix;
                    // Total acquisition time referenced to first Sample in PC memory. Demod cycle unit
                    ////////////////////////////////////////////////////////////////////////////////
                    // Case of 200% bandWidth
                    if (BandWidth == 200) {
                        TimeCycle = TimeOfFlightStartSlice + PeakDelay*DemodFrequ - PipeOffset - double(FirstSample)/4;
                        TimeCycle =TimeCycle - 0.75;            // Shift to first contributing sample for interpolation in 200% mode
                        // Test to avoidfetching data out of RF file
                        if ((TimeCycle+MinDelay) < 0) return 8;                           // ZOrigin too small for available RF data.
                        if ((TimeCycle + ((NPixSlice-1)*2*PitchPix)+MaxDelay)*4 > (NSamples-8)){
                            printf("line = %d\n", __LINE__);
                            printf("TimeCycle = %f\n", TimeCycle);
                            printf("((NPixSlice-1)*2*PitchPix) = %f\n", ((NPixSlice-1)*2*PitchPix));
                            printf("MaxDelay = %f\n", MaxDelay);
                            printf("NSamples = %d\n", NSamples);
                            printf("FirstSample = %d\n",FirstSample);
                            printf("(NSamples+FirstSample-8) = %f\n", (NSamples+FirstSample-8.0));
                            return 9;  // Image depth is too large for available RF data.
                        }
                        if (SliceStart <0) {                    // This is the first Slice contributed to by this channel
                            SliceStart = i ;                    // Index of first point contributed to by channel                   // Index of first point contributed to by channel
                            RecoLut[Indlut++] = SliceStart;                             // Indexed Lut : first RecoLut entry is start of slice
                        }
                        // Set Recolut Values according to bandwidth
                        ApodInt = int(RXApod[Arditty]*15.5);
                        if ((ApodInt >= 16) || (ApodInt < 0)) return 13;   // RXApod must be >=0 and <= 1                   // Index of first point contributed to by channel
                        if(ImmediateLut ==0) {
                            RecoLut[Indlut++] = int(32*TimeCycle+0.5);
                            RecoLut[Indlut++] = ApodInt<<3;                 // 4 higher order bits of coef index
                        }
                        else {
                            IndSample = (int(16*(TimeCycle+ReconDelay[0]*DemodFrequ)+0.5))>>2;
                            IndCoef = 2*((int(16*(TimeCycle+ReconDelay[0]*DemodFrequ)+0.5)) & 0x00000003) + 8*ApodInt;
                            RecoLut[Indlut++] = IndSample;     // Address of first data sample
                            RecoLut[Indlut++] = IndCoef;             // Immediate Lut : value of I
                        }
                    }
                    ////////////////////////////////////////////////////////////////////////////////
                    // Case of 100% bandWidth
                    else if (BandWidth == 100) {
                        TimeCycle = TimeOfFlightStartSlice + PeakDelay*DemodFrequ - PipeOffset - double(FirstSample)/2;
                        TimeCycle =TimeCycle - 1;               // Shift to first contributing sample for interpolation in 100% mode
                        // Test to avoidfetching data out of RF file
                        if ((TimeCycle+MinDelay) < 0.25) return 8;                           // ZOrigin too small for available RF data.
                        if ((TimeCycle + ((NPixSlice-1)*2*PitchPix) + MaxDelay)*2.0 > (NSamples-8.0)){
                            printf("at line = %d (TimeCycle + ((NPixSlice-1)*2*PitchPix) + MaxDelay)*2.0 > (NSamples+FirstSample-8)) \n", __LINE__);
                            printf("(TimeCycle + ((NPixSlice-1)*2*PitchPix) + MaxDelay)*2.0 = %f\n", (TimeCycle + ((NPixSlice-1)*2*PitchPix) + MaxDelay)*2.0);
                            printf("((NPixSlice-1)*2*PitchPix) = %f\n", ((NPixSlice-1)*2*PitchPix));
                            printf("MaxDelay = %f\n", (double)MaxDelay);
                            printf("NSamples = %d\n", NSamples);
                            printf("FirstSample = %d\n",FirstSample);
                            printf("transmit distance = %f\n", TransmitDistance );
                            printf("ReturnDistance = %f\n", ReturnDistance );
                            printf("XPiezo = %f\n", XPiezo );
                            printf("XSlice = %f\n", XSlice );
                            printf("ZSliceCenter = %f\n", ZSliceCenter );
                            return 9;  // Image depth is too large for available RF data.
                        }
                        if (SliceStart <0) {                    // This is the first Slice contributed to by this channel
                            SliceStart = i ;                    // Index of first point contributed to by channel
                            RecoLut[Indlut++] = SliceStart;                             // Indexed Lut : first RecoLut entry is start of slice
                        }
                        // Set Recolut Values according to bandwidth
                        ApodInt = int(RXApod[Arditty]*7.5);
                        if ((ApodInt >= 8) || (ApodInt < 0)) return 13;   // RXApod must be >=0 and <= 1
                        if(ImmediateLut ==0) {
                            RecoLut[Indlut++] = int(32*TimeCycle+0.5);
                            RecoLut[Indlut++] = ApodInt<<5;                 // 4 higher order bits of coef index
                        }
                        else {
                            IndSample = 2*((int(16*(TimeCycle+ReconDelay[0]*DemodFrequ)+0.5))>>4);
                            IndCoef = (int(16*(TimeCycle+ReconDelay[0]*DemodFrequ)+0.5)) & 0x0000000F;
                            if (IndCoef < 4 ) IndSample = IndSample-1;
                            IndCoef = 2*IndCoef + 32*ApodInt;
                            RecoLut[Indlut++] = IndSample;     // Address of first data sample
                            RecoLut[Indlut++] = IndCoef;             // Immediate Lut : value of I
                        }
                    }
                    ////////////////////////////////////////////////////////////////////////////////
                    ////////////////////////////////////////////////////////////////////////////////
                    // integrate channel contributions for central line and for 4 pixels
                    if(OutputType<100) {        //case of normalization by inverse of aperture
                        if( j== (NLines/2) ) {
                            for (p=0; p<NPixSlice; p++) {
                                AntennaGain[3*(i*NPixSlice+p)+ AntennaLeft] += float(RXApod[Arditty]);
                            }
                        }
                    }
                    else                        // Case where we want a normalization by inverse of square root of aperture
                    {
                        if( j== (NLines/2) ) {
                            for (p=0; p<NPixSlice; p++) {
                                AntennaGain[3*(i*NPixSlice+p)+ AntennaLeft] += float(RXApod[Arditty]*RXApod[Arditty]);
                            }
                        }
                    }
                }
            }
            if (SliceStart <0) {                            // No Contribution to this line
                RecoLut[Indlut++] = NSlices;                // Start beyond the end of line
            }
        }
    }

    if ( 1 == Folding ) {                // Set size of RecoLut in case of folding
        RecoLutMoinsStart[0] = Indlut;              // Offset for RecoLut on left side. In that case, size of first Lut
//         printf("RecoLutMoinsStart = %d\n", RecoLutMoinsStart[0]);
    }
    ///////////////////////////////////////////////////////
    // Compute AntennaGain array
    for (i=0; i<AntennaGainSize/2; i+=3) {
        //care of reconstructions on the right side
        if (AntennaGain[i] == 0) AntennaGain[i]= 1;                 // to avoid potential division by zero
        //care of reconstructions on the left side
        if (Folding == 0) {                                         // if no folding, specific value
            if (AntennaGain[i+1] == 0) AntennaGain[i+1]= 1;         // to avoid potential division by zero
        }
        else {             // if folding, take same value as sibling channel minus 1 if the line is centered on the reference piezo
            AntennaGain[i+1] = AntennaGain[i]- CentralChannel;
        }
//         AntennaGain[i+2]=(float)MaxGain;
//         AntennaGain[i+2]=0;
    }
    // Rearrange in SSE Form
     for (i=0; i<AntennaGainSize/2/12; i++) {
          xmma128 = _mm_set_ps(AntennaGain[12*i+ 9], AntennaGain[12*i+6], AntennaGain[12*i+3], AntennaGain[12*i  ]);
          xmmb128 = _mm_set_ps(AntennaGain[12*i+10], AntennaGain[12*i+7], AntennaGain[12*i+4], AntennaGain[12*i+1]);
          AntennaGain128[3*i  ] = xmma128;
          AntennaGain128[3*i+1] = xmmb128;
          if(OutputType<100) {           //Case of normalization by inverse of aperture
              AntennaGain128[3*i+2] = _mm_set_ps(MaxGain, MaxGain, MaxGain, MaxGain);
          }
          else {                        //Case of normalization by inverse of square root of aperture
              AntennaGain128[3*i+2] = _mm_sqrt_ps( _mm_set_ps(MaxGain, MaxGain, MaxGain, MaxGain));
          }
//           if (PlaneWave ==1) {
//              AntennaGain128[3*i+2]= _mm_min_ps(AntennaGain128[3*i+2], _mm_sqrt_ps( _mm_set_ps(A+(4*i+3)*B, A+(4*i+2)*B, A+(4*i+1)*B, A+(4*i+0)*B)));
//           }
    }
    ////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////////////////////////
    // Computing globalLut
    for (NRec = 0; NRec < NRecon; NRec++) {        //loop on all reconstructions
        if(ReconAbsc[NRec] < ReconMux[NRec]) return 15;  // Region cannot be left of aperture
        if(ReconAbsc[NRec] >= (ReconMux[NRec] + NChannels)) return 16;  // Region cannot be right of aperture
        GlobalLut[NRec] = (ReconMux[NRec]& 0x000000FF) + 256*((ReconAbsc[NRec]) & 0x000000FF) + 256*256*int(32*ReconDelay[NRec]*DemodFrequ + 0.5);
    }
//     printf("RecoLutMoinsStart = %d\n", RecoLutMoinsStart[0]);

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Transmit 'AcqInfo' and 'ReconInfo' structures
    // 'Probe'
    ProbeOut[0] = NPiezos;              // integer unit. number of piezos on probe
    ProbeOut[1] = PiezoWidth;           // mm unit. Piezo width.
    ProbeOut[2] = PiezoPitch;           // mm unit. Piezo pitch.
    // 'AcqInfo'
    AcqInfoOut[0] = SteerAngle;          // radian unit. Steering angle.
    AcqInfoOut[1] = BandWidth;            // 200 100 50. (I,Q) acquisition bandwidth.
        // FrequDivider not used
    AcqInfoOut[3] = DemodFrequ;          // MHz unit. Demodulation frequency.
    AcqInfoOut[4] = SyntAcq;             // Synthetic acquisition flag. 0 or 1.
    AcqInfoOut[5] = ZOrigin;             // mm unit. depth of first reconstructed points.
        // Depth not used
        // InitialTime not used
    AcqInfoOut[8] = FirstSample;          // integer. Index of first sample sent to PC (useful for partial transmission).
    AcqInfoOut[9] = NSamples;             // integer. Number of samples per channel.
        // NAcqChannels not used
    // 'ReconInfo'
    ReconInfoOut[ 0] = NRecon;               // integer unit. Number of reconstructions
    ReconInfoOut[ 1] = PlaneWave;            // 0 or 1. 1 if plane wave reconstruction, 0 if focalized reconstruction
    ReconInfoOut[ 2] = NLines;              // integer unit. Number of lines in region. 1 for plane wave, >= 1 for focused wave.
    ReconInfoOut[ 3] = XOrigin;             // PiezoPitch unit. Left line abscissa offset/left edge of reference piezo. usually 0 for plane Wave.
    ReconInfoOut[ 4] = LinePitch;           // PiezoPitch unit. Horizontal distance between two lines. 1 for normal resolution, 0.5 for double resolution.
    ReconInfoOut[ 5] = PitchPix;            // Lambda unit. Distance between two pixels on reconstructed line.
    ReconInfoOut[ 6] = NPixLine;             // integer unit. Number of reconstructed slices
    ReconInfoOut[ 7] = PeakDelay;            // µsec unit . Delay between start of transmit pulse, and maximum of transmit power.
    ReconInfoOut[ 8] = PipeOffset;           // integer. demod cycle unit . Time corresponding to first sample, assumed to be a 'I'.
    ReconInfoOut[ 9] = NChannels;           // integer. Number of acquisition channels.
    ReconInfoOut[10] = MaxChannelsOneSide;  // Maximum number of channels contributing on one side of a region.
    ReconInfoOut[11] = ChannelOffset;       // Channel offset in receive buffer.
    ReconInfoOut[12] = TxAlongRx;           // 1 = Receive aperture is set around firing abscissa.
    ReconInfoOut[13] = BandInterp;           // real [0 100]. Low pass filter for interpolation. 97 or 65
    ReconInfoOut[14] = FootInterp;           // real . 100 * Foot of Hanning interpolation. 50 or 14
    ReconInfoOut[15] = RXApodOne;           // max aperture value for full gain (Arditty coefficient)
    ReconInfoOut[16] = RXApodZero;          // max aperture value for a contributing channel
    ReconInfoOut[17] = OutputType;          // Type of output : 0=(I,Q)   1= intensity
    ReconInfoOut[18] = NPixSlice;           // integer unit. Number of reconstructed pixels per slice.
    ReconInfoOut[19] = LutSizeBytes;        // Lut size in bytes
    ReconInfoOut[20] = Mode;                 // 0 or 1. 1 if folding is to be performed, else 0.
    ReconInfoOut[21] = SoundSpeed;          // m/sec unit. Speed of sound.

  //////////////////////////////////////////////////////////////////////////////////////////////////
    return 0;

} //end of function
////////////////////////////////////////////////////////////////////////////////////////////////////
