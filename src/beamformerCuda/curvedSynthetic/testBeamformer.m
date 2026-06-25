clear bfStruct Data img

bfStruct.fNumber = 0.8 ; % fnumber
bfStruct.xOrigin = 0.5;
bfStruct.nbLinesPerRecon = 1;

bfStruct.ProbeRadius =system.probe.Radius*1e-3 ;
bfStruct.speedOfSound = common.constants.SoundSpeed;
bfStruct.piezoPitch = system.probe.Pitch*1e-3;
bfStruct.demodFreq = UF.TwFreq*1e6;
bfStruct.rOrigin = 5*1e-3;
bfStruct.peakDelay = 0;
bfStruct.nbAngles = length(UF.FlatAngles);
bfStruct.reconPitch = UFBfInfo.LinePitch*system.probe.Pitch*1e-3;
bfStruct.nbRecon = UF.RxWidth*1e-3/bfStruct.reconPitch;

bfStruct.linePitch = bfStruct.reconPitch ;
bfStruct.lambda = bfStruct.speedOfSound/bfStruct.demodFreq;
bfStruct.pixelPitch = UFBfInfo.PitchPix*bfStruct.lambda;
bfStruct.nbPiezos = system.probe.NbElemts;
bfStruct.firstSample = Sequence.InfoStruct.event(1).skipSamples;
bfStruct.nbSamples = Sequence.InfoStruct.event(1).numSamples;
bfStruct.nbPixelsPerLine = (UFBfInfo(1).Depth(2)-UFBfInfo(1).Depth(1))/bfStruct.pixelPitch*1e-3;
bfStruct.normMode=2;
bfStruct.synthAcq = (UF.TxWidth>128*system.probe.Pitch);
bfStruct.usegpu = 1;
bfStruct.firstAcquiredChannel =  max(0,min(system.probe.NbElemts,round(UF.RxCenter- UF.RxWidth/2)/system.probe.Pitch));
bfStruct.firstRecon = (UF.TxCenter - UF.TxWidth/2)/system.probe.Pitch;
bfStruct.angles = UF.FlatAngles*pi/180 ;
bfStruct.timeOrigin = abs(bfStruct.nbPiezos*bfStruct.piezoPitch/2*sin(bfStruct.angles)/bfStruct.speedOfSound);
bfStruct.channelOffset = Sequence.InfoStruct.mode(1).channelSize/2;

tic
IQ = beamformerCurvedSynthetic(bfStruct,buffer.data,buffer.alignedOffset);
toc
img = abs(complex(IQ(1:2:end,:),IQ(2:2:end,:)));
imagesc(20*log10(abs(img)));

angProbePitch = system.probe.Pitch/system.probe.Radius ;
lambda = common.constants.SoundSpeed / bfStruct.demodFreq ;

parameters.NR = size(img,1);
parameters.NTheta = size(img,2);
parameters.Nx = 512 ;
parameters.Nz = 512 ;

% TODO: FirstLinePos is in number of elements and XOrigin is in mm
parameters.AngleOfFirstLine = (ImgInfo.FirstLinePos + ImgInfo.XOrigin - 1) * angProbePitch - system.probe.NbElemts * angProbePitch/2 ;
parameters.AngleStep = bfStruct.reconPitch*1e3/system.probe.Radius ;
parameters.FirstSampleAxialPosition = system.probe.Radius + bfStruct.rOrigin*1e3 ;
parameters.RadialStep = bfStruct.pixelPitch * 1e3 ;
parameters.OffSet = zeros(2,1, 'double');
parameters.SizeOfImg = zeros(2,1, 'double');
parameters.NbVois = 1 ;
parameters.SmoothDistance = 1e3*lambda/2;
parameters.SecondPass = 0 ;
parameters.OffSetX = 0 ;
parameters.OffSetZ = 0 ;
parameters.SizeOfImgX = 0 ;
parameters.SizeOfImgZ = 0 ;
parameters.ComputeOffSetAndSize = 1 ;

lutsize = zeros(1,'int32');

parameters_array = processing.SCV_Parser( parameters ) ;
processing.ScanConversionTable( parameters_array, lutsize );

parameters.OffSetX = parameters_array(10);
parameters.OffSetZ = parameters_array(11);
parameters.SizeOfImgX = parameters_array(12);
parameters.SizeOfImgZ = parameters_array(13);

NbVoistot = ( parameters.NbVois * 2 ) ^ 2 ;

parameters.SecondPass = 1 ;
LUT_coeff = zeros( lutsize*NbVoistot, 1, 'single' );
LUT_X     = zeros( lutsize, 1, 'int16' );
LUT_Z     = zeros( lutsize, 1, 'int16' );
LUT_Vois1 = zeros( lutsize*NbVoistot, 1, 'int16' );
LUT_Vois2 = zeros( lutsize*NbVoistot, 1, 'int16' );

processing.ScanConversionTable( processing.SCV_Parser(parameters), LUT_coeff, LUT_X, LUT_Z, LUT_Vois1, LUT_Vois2, lutsize );

Data.XAxis = parameters.OffSetX + (0:(parameters.Nx-1)) * parameters.SizeOfImgX / (parameters.Nx-1) ;
Data.ZAxis = parameters.OffSetZ + (0:(parameters.Nz-1)) * parameters.SizeOfImgZ / (parameters.Nz-1) ;
Data.OUT   = zeros(parameters.Nx,parameters.Nz,'single');
Data.LUT_coeff = LUT_coeff ;
Data.LUT_X = LUT_X ;
Data.LUT_Z = LUT_Z ;
Data.LUT_Vois1 = LUT_Vois1 ;
Data.LUT_Vois2 = LUT_Vois2 ;
Data.parameters = parameters;
Data.lutsize = lutsize ;
Data.img = zeros(Nz,Nx);

Data = processing.ScanConvertImage( Data, img );
X = Data.XAxis;
Z = Data.ZAxis;
%%
figure(2);
imagesc( X, Z, Data.OUT', [70 120] )
colormap(gray); colorbar;
axis equal; axis tight;
t = toc ;
drawnow ;

