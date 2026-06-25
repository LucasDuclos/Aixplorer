#include "mex.h"
#include "reconstruct.h"

#if defined _WIN32 || defined _WIN64
    DWORD WINAPI spWinThreadFunction( LPVOID lpParam ) 
    {
        return (DWORD)sp::reconstruct ((void *)lpParam);
    }
#endif

/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
    bool isComplexIQData ;
    int NThreads, thrd ;
    int Status ;
    mxArray * LutStruct ;
    mxArray * buffer ; 
    mxArray * IQData ;
    mxArray * RFData ;
    short *RcvDataShort ;
    float * IQPix ;
    float * IQPixImag = NULL ;
    mxArray * pArray ;
      
    int DebugMode ;
    
    // Check for proper number of arguments.
    if (nrhs != 4)
        mexErrMsgTxt("Four inputs required.");
    if (nlhs != 0)
        mexErrMsgTxt("Zeros output required.");
    
    
    // get inputs
    LutStruct =(mxArray*)(prhs[0]) ; // get lut structures
    NThreads =  mxGetNumberOfElements(LutStruct); // NThread = NLuts
    
    if(NThreads==0)
    {
        mexErrMsgTxt("empty lut structure");
    }
    
    buffer = (mxArray*)(prhs[1]) ; // get buffer
    
    if(mxGetNumberOfElements(buffer)!=1)
    {
        mexErrMsgTxt("buffer structure must be of length 1");
    }
    // get buffer data
    
    pArray = mxGetField(buffer, (int)(0), "data") ;
    
    if(pArray && !mxIsEmpty(pArray))
    {
        
        RFData = (mxArray*)pArray ;
    }
    else
    {
        mexErrMsgTxt("buffer structure does not contains requiered 'data' field or this field is empty\n");
    }
    
    // get alignedOffset
    pArray = mxGetField(buffer, (int)(0), "alignedOffset") ;
    int alignedOffset = 0 ;
    if(pArray && !mxIsEmpty(pArray))
    {
        alignedOffset = (int)(*(unsigned int*)mxGetData(pArray))/2 ;
        //mexPrintf("alignedOffset = %d\n",alignedOffset);
    }
    else
    {
        mexErrMsgTxt("buffer structure does not contains requiered 'alignedOffset' field or this field is empty\n");
    }
           
    IQData = (mxArray*)(prhs[2]); // get IQ Data output matrix
    DebugMode = (int)mxGetScalar(prhs[3]);
    
    // get info on RFbuffer
    if(mxGetClassID(RFData)!=mxINT16_CLASS)
    {
        mexErrMsgTxt("RF data must be int16");
    }
    
    int totalNbRFSamples = mxGetNumberOfElements(RFData);
    RcvDataShort = (short*)mxGetData(RFData);
    
    // get info and check output buffer    
    if(mxGetClassID(IQData)!=mxSINGLE_CLASS)
    {
        mexErrMsgTxt("IQ data must be single");
    }
    
    isComplexIQData = mxIsComplex(IQData) ;
    if(isComplexIQData)
    {
        IQPix = (float*)mxGetData(IQData);
        IQPixImag = (float*)mxGetImagData(IQData);
    }
    else
    {
        IQPix = (float*)mxGetData(IQData);
        IQPixImag = NULL ;
    }
    
    // get info from IQData buffer
    int nbDims = mxGetNumberOfDimensions(IQData);
    int * Dims = new int[nbDims];
    int Nz, Nx, Ni, Npixels, NpixelsPerImage ;
    
    if(nbDims==2) {
        Nz = mxGetDimensions(IQData)[0];
        Nx = mxGetDimensions(IQData)[1];
        Ni = 1 ;
    }
    else if(nbDims==3) {
        Nz = mxGetDimensions(IQData)[0];
        Nx = mxGetDimensions(IQData)[1];
        Ni = mxGetDimensions(IQData)[2];
    }
    else {
        mexErrMsgTxt("IQData must have 2 or 3 dimensions");
    }
    
    Npixels = Nz*Nx*Ni ; // total number of pixels    
    NpixelsPerImage = Nz*Nx ; // number of pixels per image
    
        
    if(DebugMode>0) // Debug test
    {
        mexPrintf("DebugMode = %d\n",DebugMode);
        mexPrintf("Nz = %d, Nx = %d, Ni = %d \n",Nz,Nx,Ni);        
        int ProbeSize = 3*4;
        int AcqInfoSize = 11*4;
        double * ReconInfo, * AcqInfo ;
        
        // tests to avoid leaks
        
        //int * firstPixel = new int [NThreads];
        //int * lastPixel = new int [NThreads];
        
        int   NPixLine = 0;
        int   NLines = 0;
        int   NSamples = 0;
        int   SyntAcq = -1 ;
        
        int * NRecon = new int[NThreads];
        int * FrameID = new int[NThreads];
        int * LineID = new int[NThreads];
        int * OutputType = new int[NThreads];
        int * firstPixel = new int[NThreads];
        int * lastPixel = new int[NThreads];
        int * PlaneWave = new int[NThreads];
        
        bool isFrameIDgt0 = false ;
        bool isLineIDgt0 = false ;
       
        
        for(int lutIdx=0;lutIdx<NThreads;lutIdx++)
        {
                        
            // FrameID
            pArray = mxGetField(LutStruct, lutIdx, "FrameId") ;
            if(pArray && !mxIsEmpty(pArray))
            {
                FrameID[lutIdx] = (int)(mxGetScalar(pArray));
                if(FrameID[lutIdx]>0)
                {
                   isFrameIDgt0 = true ; 
                }
            }
            else
            {
                mexPrintf("lut Id = %d\n",lutIdx);
                mexErrMsgTxt("FrameId field does not exist or is empty in lut");
            }
            
            // LineID
            pArray = mxGetField(LutStruct, lutIdx, "LineId") ;
            if(pArray && !mxIsEmpty(pArray))
            {
                LineID[lutIdx] = (int)(mxGetScalar(pArray)) ;
                if(LineID[lutIdx]>0)
                {
                   isLineIDgt0 = true ; 
                }
            }
            else
            {
                mexPrintf("lut Id = %d\n",lutIdx);
                mexErrMsgTxt("LineID field does not exist or is empty in lut");
            }
            
            // test of Line and Frame IDs
            if(LineID[lutIdx]>0 && FrameID[lutIdx]>0)
            {
                mexPrintf("lut Id = %d\n",lutIdx);
                mexErrMsgTxt("LineID and FrameID cannot be >0 at the same time, error in lut");
            }
            
            // Lut
            short * Lut = NULL ;
            pArray = mxGetField(LutStruct, lutIdx, "data") ;
            if(pArray && !mxIsEmpty(pArray))
            {
                Lut = (short*)mxGetData(pArray);
            }
            else
            {
                mexPrintf("lut Id = %d\n",lutIdx);
                mexErrMsgTxt("data field does not exist or is empty in lut");
            }
                       
            // get lut info
            ReconInfo = ( double * ) &Lut[ProbeSize + AcqInfoSize];
            AcqInfo = ( double * ) &Lut[ProbeSize];

            // parse and control NSamples
            if(NSamples==0)
            {
                NSamples = ( int ) AcqInfo[9] ;
            }
            else
            {
                if(NSamples != ( int ) AcqInfo[9])
                {
                    mexPrintf("lut Id = %d\n",lutIdx);
                    mexErrMsgTxt("NSamples must be identical for all luts");
                }
            }
            
            // parse and control SyntAcq
            if(SyntAcq==-1)
            {
                SyntAcq = ( int ) AcqInfo[4];
            }
            else
            {
                if(SyntAcq != ( int ) AcqInfo[4])
                {
                    mexPrintf("lut Id = %d\n",lutIdx);
                    mexErrMsgTxt("SyntAcq must be identical for all luts");
                }
            }            
            
            // parse and control NPixLine
            if(NPixLine==0)
            {
                NPixLine = ( int ) ReconInfo[6] ;  // integer unit. Number of reconstructed slices
            }
            else
            {
                if(NPixLine != ( int ) ReconInfo[6] )
                {
                    mexPrintf("lut Id = %d\n",lutIdx);
                    mexErrMsgTxt("NPixLine must be identical for all luts");
                }
            }
            
            // parse and control NPixLine
            if(NLines==0)
            {
                NLines = ( int ) ReconInfo[2] ;
            }
            else
            {
                if(NLines != ( int ) ReconInfo[2] )
                {
                    mexPrintf("lut Id = %d\n",lutIdx);
                    mexErrMsgTxt("NLines must be identical for all luts");
                }
            }
            
            
            // variable params
            NRecon[lutIdx] =  (int)ReconInfo[0];                                 
            OutputType[lutIdx] = (int)ReconInfo[17];
            PlaneWave[lutIdx] = ( int ) ReconInfo[1] ;

            // first and last pixel writed in IQData for this lut (used later tests)
            firstPixel[lutIdx] = LineID[lutIdx]*NPixLine + FrameID[lutIdx]*NPixLine*NLines*NRecon[lutIdx] ;
            lastPixel[lutIdx] = LineID[lutIdx]*NPixLine + (FrameID[lutIdx]+1)*NPixLine*NLines*NRecon[lutIdx] -1 ;
            
            // a stupid test
            if(lastPixel[lutIdx]<firstPixel[lutIdx])
            {
                mexPrintf("lut Id = %d\n",lutIdx);
                mexErrMsgTxt("zero pixel in lut");
            }            
            
        } // end for (
        
        // test 1: FrameID/LineID
        
        if(isLineIDgt0 && isFrameIDgt0)
        {
            mexErrMsgTxt("FrameID and LineID >0");
            
           // if( (NPixLine != Nz) || ( (NLines*NReconAll + LineID) != Nx) )
           //         mexErrMsgTxt("IQ data matrix is to small");
            
        }
        else if(isLineIDgt0 && !isFrameIDgt0) // bmode, all the thread work on the same frame
        {
            int NReconAll = 0;
            for(int lutIdx=0;lutIdx<NThreads;lutIdx++)
            {
                NReconAll+=NRecon[lutIdx] ;

                // check individually allocation
                if( NPixLine*NLines*NReconAll > Npixels )
                    mexErrMsgTxt("IQ data matrix is to small");
                
                // check complexity
                if(!isComplexIQData && ( (OutputType[lutIdx]%100) == 3 && (OutputType[lutIdx]%100) == 0))
                    mexErrMsgTxt("IQ data must be complex because (I,Q) output is requested as OutputType mod(10) = 0 or 3 ");
                
                if(OutputType[lutIdx]==13)
                {
                        mexErrMsgTxt("Forbidden output type (13)");
                }
                    
                // check for potential overlap when writing pixels
                for(int lutIdx2=lutIdx+1;lutIdx2<NThreads;lutIdx2++)
                {
                    bool test1 =  (firstPixel[lutIdx]>= firstPixel[lutIdx2]) && (firstPixel[lutIdx] <= lastPixel[lutIdx2])  ;
                    bool test2 =  (lastPixel[lutIdx]>= firstPixel[lutIdx2]) && (lastPixel[lutIdx] <= lastPixel[lutIdx2])  ;
                    
                    if(test1 || test2)
                    {
                        mexPrintf("lut Id = %d and %d\n",lutIdx,lutIdx2 );
                        mexPrintf("firstPixel[%d] = %d,lastPixel[%d]= %d\n", lutIdx,firstPixel[lutIdx],lutIdx,lastPixel[lutIdx]);
                        mexPrintf("firstPixel[%d] = %d,lastPixel[%d]= %d\n", lutIdx2,firstPixel[lutIdx2],lutIdx2,lastPixel[lutIdx2]);
                        mexErrMsgTxt(" pixels writed are overlaping");
                    }
                }
            }
            // check parameters are consistant
                if( (NPixLine != Nz) || ( (NLines*NReconAll ) != Nx) )
                    mexErrMsgTxt("IQ data matrix is to small, check dimension 1 and 2 are consistant with NPixLine and NLines, NRecon");
        }
        else // FrameID >=0 and Lineid = 0
        {
            int NReconAll = (int)NRecon[0];
            for(int lutIdx=0;lutIdx<NThreads;lutIdx++)
            {

                // check NRecon
                if(NReconAll != NRecon[lutIdx])
                {
                    mexPrintf("lut Id = %d\n",lutIdx);
                    mexErrMsgTxt("NRecon must be the same for each Lut");
                }
                // check individual allocation
                if( NPixLine*NLines*NReconAll*(FrameID[lutIdx]+1) > Npixels )
                    mexErrMsgTxt("IQ data matrix is to small");

                // check parameters are consistant
                if( (NPixLine != Nz) || ( (NLines*NReconAll ) != Nx) )
                    mexErrMsgTxt("IQ data matrix is to small, check dimension 1 and 2 are consistant with NPixLine and NLines, NRecon");

                // check complexity
                if(!isComplexIQData && ( (OutputType[lutIdx]%100) == 3 && (OutputType[lutIdx]%100) == 0))
                    mexErrMsgTxt("IQ data must be complex because (I,Q) output is requested as OutputType mod(10) = 0 or 3 ");
                
                if(OutputType[lutIdx]==13)
                {
                        mxArray * pTemp = mxGetField(LutStruct, thrd , "IQbuffer");
                        if( (pTemp == NULL ) || ( mxGetNumberOfElements(pTemp) < (2*NPixLine*NReconAll*NLines)) )
                        {
                            mexErrMsgTxt("IQbuffer is too small or does not exist");
                        }
                }
                    
                // check for potential overlap when writing pixels
                for(int lutIdx2=lutIdx+1;lutIdx2<NThreads;lutIdx2++)
                {
                    bool test1 =  (firstPixel[lutIdx]>= firstPixel[lutIdx2]) && (firstPixel[lutIdx] <= lastPixel[lutIdx2])  ;
                    bool test2 =  (lastPixel[lutIdx]>= firstPixel[lutIdx2]) && (lastPixel[lutIdx] <= lastPixel[lutIdx2])  ;
                    
                    if(test1 || test2)
                    {
                        mexPrintf("lut Id = %d and %d\n",lutIdx,lutIdx2 );
                        mexErrMsgTxt(" pixels writed are overlaping");
                    }
                }
            }
           
        } //else // FrameID >=0 and Lineid = 0
        
        if(DebugMode>1) // check RF buffer
        {
            
            int *  firstRFSampleperChannel = new int[NThreads];
            int *  lastRFSampleperChannel = new int[NThreads];
            
            for(int lutIdx=0;lutIdx<NThreads;lutIdx++)
            {                
                mxArray * pTemp = NULL ;
                
                pTemp = mxGetField(LutStruct, lutIdx, "RcvOffset") ;
                
                firstRFSampleperChannel[lutIdx] = 0 ;
                if(pTemp) {
                    firstRFSampleperChannel[lutIdx] = (int)mxGetScalar(pTemp);
                }
                else {
                    mexPrintf("lut Id = %d\n", lutIdx);
                    mexErrMsgTxt("RcvOffset field missing");
                }
                
                if(PlaneWave == 0) // focused mode, skip one RF data frame per Recon (xNRecon)
                    lastRFSampleperChannel[lutIdx] = firstRFSampleperChannel[lutIdx] + (1+SyntAcq)*NRecon[lutIdx]*NSamples -1 ;
                else // plane wave of other (PAT must be ok but it need check for RF accumulation !!!)
                    lastRFSampleperChannel[lutIdx] = firstRFSampleperChannel[lutIdx] + (1+SyntAcq)*NSamples -1 ;
                
            }
            
            // check for overlap !
            for(int lutIdx=0;lutIdx<NThreads;lutIdx++)
            {                
                // check if RF buffer is big enough
                
                if(lastRFSampleperChannel[lutIdx]>=totalNbRFSamples)
                {
                    mexPrintf("lut Id = %d\n",lutIdx );
                    mexErrMsgTxt(" stop here to avoid segmentation fault : idx RF sample > length(RFData(:))");
                }                    
                
                // check for potential overlap when writing pixels
                for(int lutIdx2=lutIdx+1;lutIdx2<NThreads;lutIdx2++)
                {
                    bool test1 =  (firstRFSampleperChannel[lutIdx]>= firstRFSampleperChannel[lutIdx2]) && (firstRFSampleperChannel[lutIdx] <= lastRFSampleperChannel[lutIdx2])  ;
                    bool test2 =  (lastRFSampleperChannel[lutIdx]>= firstRFSampleperChannel[lutIdx2]) && (lastRFSampleperChannel[lutIdx] <= lastRFSampleperChannel[lutIdx2])  ;
                    
                    if(test1 || test2)
                    {
                        mexPrintf("lut Id = %d and %d",lutIdx,lutIdx2 );
                        mexErrMsgTxt(" RF sample readed are overlaping");
                    }
                }
                
            }
            
            delete [] firstRFSampleperChannel;
            delete [] lastRFSampleperChannel;
            
        } // end if debug mode > 1
                      
        delete [] NRecon ;
        delete [] FrameID ;
        delete [] LineID ;
        delete [] OutputType ;
        delete [] firstPixel ;
        delete [] lastPixel ;
        delete [] PlaneWave ;
                 
	} // end if debug
        
    ArgReconstruct *Arg = new ArgReconstruct[NThreads];
    
    
    for(thrd=0;thrd<NThreads;thrd++) {
        
        
        // offset RF
        
        mxArray * pTemp = NULL ;
        
        pTemp = mxGetField(LutStruct, thrd, "RcvOffset") ;
        int inputoffset = 0 ;
        if(pTemp)
        {
            inputoffset = (int)mxGetScalar(pTemp) + alignedOffset ;
        }
        else
        {
             mexPrintf("lut Id = %d ",thrd);
             mexErrMsgTxt("RcvOffset field missing");            
        }

        // offset on the first pixel of the image to be beamformed
        pTemp = mxGetField(LutStruct, thrd, "FrameId") ;
        int outputOffset = 0 ;
        if(pTemp)
        {
             outputOffset = NpixelsPerImage*(int)mxGetScalar(pTemp) ; // frameId
        }
        else
        {
             mexPrintf("lut Id = %d ",thrd);
             mexErrMsgTxt("FrameId field missing");            
        }
        
        pTemp = mxGetField(LutStruct, thrd, "LineId") ;
        if(pTemp)
        {
             outputOffset += Nz*(int)mxGetScalar(pTemp) ; //lineId
        }
        else
        {
             mexPrintf("lut Id = %d ",thrd);
             mexErrMsgTxt("LineId field missing");            
        }                

        // lut data
        pTemp = mxGetField(LutStruct, thrd, "data");
        if(pTemp)
        {
             Arg[thrd].Lut = (short*)mxGetData(pTemp);
        }
        else
        {
             mexPrintf("lut Id = %d ",thrd);
             mexErrMsgTxt("data field missing");            
        }               
        
        // IQ data output buffers
        Arg[thrd].IQbuffer = NULL;
        pTemp = mxGetField(LutStruct, thrd , "IQbuffer");
        if(pTemp)
        {
            Arg[thrd].IQbuffer = (float*)mxGetData(pTemp);            
        }
        else
        {
            Arg[thrd].IQbuffer = NULL;
        }
        
        // final test on RF buffer size
        
        double * arrayParam = (double*)(Arg[thrd].Lut);
        int ProbeSize = 3;
        int AcqInfoSize = 11;
        int ReconInfoSize = 26;
        int NSamples = (int) (arrayParam[ProbeSize+9]) ;
        int channelSize = ( int ) ( arrayParam[ProbeSize + AcqInfoSize+11] /2) ;
        int SyntAcq = ( int ) arrayParam[ProbeSize +4]; 
        int RxChan = (int) arrayParam[ProbeSize + AcqInfoSize+9]; // taking into account syntAcq
        if((inputoffset - alignedOffset +  (1+SyntAcq)*NSamples ) > channelSize)
        {
            mexPrintf("RcvOffset = %d, alignedOffset = %d, NSamples = %d, SyntAcq = %d, channelSize = %d\n",
                 inputoffset-alignedOffset,alignedOffset,NSamples,SyntAcq,channelSize);
            mexErrMsgTxt("aborting reconstruction to prevent segmentation, check RF offset");
        }
        
        int nbRFSamplesAfterOffset = totalNbRFSamples - inputoffset ;
        int nbExpectedRFSamples = 0 ;
        int PlaneWave = (int) arrayParam[ProbeSize + AcqInfoSize + 1];
        if(PlaneWave ==0)
        {
            int NRecon = (int) arrayParam[ProbeSize + AcqInfoSize] ;
            nbExpectedRFSamples = NRecon*NSamples*RxChan ;
        }
        else
        {
            nbExpectedRFSamples = NSamples*RxChan ;
        }
        
        if(nbRFSamplesAfterOffset<nbExpectedRFSamples)
        {
            mexPrintf("nbRFSamplesAfterOffset = %d,nbExpectedRFSamples = %d\n",  nbRFSamplesAfterOffset,nbExpectedRFSamples);
            mexErrMsgTxt("aborting reconstruction to prevent segmentation, check RF offset");            
        }
                
        // input and output buffer pointers
        Arg[thrd].IQPix = IQPix + outputOffset ;
        Arg[thrd].RcvDataShort = &(RcvDataShort[inputoffset]);
        Arg[thrd].dataFormat = E_MATLAB ;
        if(isComplexIQData)
        {
            Arg[thrd].IQPixImag = IQPixImag +outputOffset ;
        }
    }
    
    #ifdef __linux__
        pthread_t *handle = new pthread_t[NThreads] ;

        // launch thread
        for(thrd=0;thrd<NThreads;thrd++) {
            pthread_create(&handle[thrd], NULL, sp::reconstruct, &Arg[thrd]);
        }

        for(thrd=0;thrd<NThreads;thrd++) {
            pthread_join(handle[thrd], NULL);
        }
        delete [] handle;
     
    #elif defined _WIN32 || defined _WIN64
        HANDLE  *hThreadArray = new HANDLE[NThreads];
        DWORD   *dwThreadIdArray = new DWORD[NThreads];
        for(thrd=0;thrd<NThreads;thrd++)
        {
            hThreadArray[thrd] = CreateThread( NULL,                   // default security attributes
                                               0,                      // use default stack size  
                                               spWinThreadFunction,        // thread function name
                                               &Arg[thrd],             // argument to thread function 
                                               0,                      // use default creation flags 
                                               &dwThreadIdArray[thrd]);   // returns the thread identifier  
        }
        WaitForMultipleObjects(NThreads,hThreadArray,TRUE,INFINITE);
        delete [] hThreadArray;
        delete [] dwThreadIdArray;
    #else
        #error "unknown platform"
    #endif
    

    delete [] Arg ;
}
