% SCRIPTS.ULTRAFASTCURVEDAPEX (PUBLIC)
%   Build and run a sequence with an ULTRAFASTCURVEDAPEX imaging.
%
%   Note - This function is defined as a script of SCRIPTS package. It cannot be
%   used without the legHAL package developed by SuperSonic Imagine and without
%   a system with a REMOTE server running.
%
%   Copyright 2016 Supersonic Imagine
%   Revision: 1.00 - Date: 2010/03/25

%% Parameters definition

% System parameters
debugMode = 0;
AixplorerIP    = '127.0.0.1';% IP address of the Aixplorer device
ImagingVoltage = 50;            % imaging voltage [V]
ImagingCurrent = 1;             % security current limit [A]

% ============================================================================ %

repetition = 1 ;

% Image parameters
ImgInfo.Depth(1)           = 2;% initial depth [mm]
ImgInfo.Depth(2)           = 120; % image depth [mm]
ImgInfo.CorrectFactorDepth = 0.2;% depth correction to estimate the acquisition duration
ImgInfo.NbColumnsPerPiezo  = 1;
ImgInfo.NbRowsPerLambda    = 1;

% =======================================================================²==== %

% ULTRAFASTCURVEDAPEX acquisition parameters

angulator = 0;

%% Compute apex positions

% SC6-1 Parameters for libssimath utest comparison
% system.probe.Pitch = 0.33;
% system.probe.Radius = 60;
% system.probe.NbElemts = 192;
box = tools.polarbox(1.5*system.probe.Radius, pi/12, ...
                        0.75*system.probe.Radius, pi/9);
                                     
% deltaThetaProbe =  pi/6;
deltaThetaProbe =  system.probe.Pitch*system.probe.NbElemts/(2.0*system.probe.Radius);
probe = tools.polarprobe(system.probe.Radius, deltaThetaProbe);

UF.thApex       = (-9:3:9+angulator)*pi/180.;

for i = 1:1:length(UF.thApex);
    
    apx = tools.apex(box, probe, UF.thApex(i));

    UF.xApex(i) = apx.xI; 
    UF.zApex(i) = apx.zI;
    
    [TxCenter(i), TxWidth(i)] = apx.computeTXAperture(box, probe);
    
end


UF.PulseInversion = 0;            % 1 or 0
UF.TxFreq       = 45/16 ;         % emission frequency [MHz]
UF.NbHcycle     = 4;              % # of half cycles
% UF.xApex        = -UF.rApex*sin(UF.thApex);
% UF.zApex        = -UF.rApex*cos(UF.thApex)+system.probe.Radius;
UF.DutyCycle    = 1;              % duty cycle [0, 1]
UF.ApodFct      = 'none';         % emission apodization function (none, hanning)
UF.RxBandwidth  = 2;              % sampling mode (1 = 200%, 2 = 100%, 3 = 50%)
UF.FIRBandwidth = 40;             % FIR receiving bandwidth [%] - center frequency = UFTxFreq
UF.NbFrames     = 1;              % # of acquired images
UF.PRF          = 400;            % pulse frequency repetition [Hz] (0 for greatest possible)
UF.FrameRate    = 0;              % pulse frequency repetition [Hz] (0 for greatest possible)
UF.TGC          = [900 900 900 900 900 900]; % TGC profile
UF.TrigIn       = 0;

% ============================================================================ %
% ============================================================================ %

%% DO NOT CHANGE - Additional parameters

% Additional parameters for ULTRAFASTPHASED acquisition
UF.RxFreq = (1+UF.PulseInversion) * 4 * UF.TxFreq;% sampling frequency [MHz]

% Estimate RxDelay for ULTRAFASTPHASED acquisition
UF.RxDelay = floor((1 - ImgInfo.CorrectFactorDepth) * ImgInfo.Depth(1) * 1e3 ...
    * UF.RxFreq / common.constants.SoundSpeed) / UF.RxFreq;

% Estimate RxDuration for ULTRAFASTPHASED acquisition
UF.RxDuration = ceil((1 + ImgInfo.CorrectFactorDepth) ...
    * (2 * ImgInfo.Depth(2) - ImgInfo.Depth(1)) * 1e3 ...
    * (UF.RxFreq / 2^(UF.RxBandwidth - 1)) ...
    / (common.constants.SoundSpeed * 128)) * 128 ...
    / (UF.RxFreq / 2^(UF.RxBandwidth - 1)) ;

%UF.RxDuration = UF.RxDuration*2;

% ============================================================================ %
% ============================================================================ %

%% DO NOT CHANGE -  acquisition mode

try
    
    % Create the ULTRAFASTPHASED acquisition mode and add the push elusev
    Acmo = acmo.ultraphased( ...
        'TwFreq',       UF.TxFreq, ...
        'NbHcycle',     UF.NbHcycle, ...
        'xApex',        UF.xApex, ...
        'zApex',        UF.zApex, ...
        'DutyCycle',    UF.DutyCycle, ...
        'TxCenter',     TxCenter, ...
        'TxWidth',      TxWidth, ...
        'ApodFct',      UF.ApodFct, ...
        'RxFreq',       UF.RxFreq, ...
        'RxDuration',   UF.RxDuration, ...
        'RxDelay',      UF.RxDelay, ...
        'RxCenter',     system.probe.NbElemts*system.probe.Pitch/2 + angulator*system.probe.Radius/180*pi, ...
        'RxWidth',      system.probe.NbElemts*system.probe.Pitch, ...
        'RxBandwidth',  UF.RxBandwidth, ...
        'FIRBandwidth', UF.FIRBandwidth, ...
        'RepeatFlat',   UF.NbFrames, ...
        'PRF',          UF.PRF, ...
        'FrameRate',    UF.FrameRate, ...
        'ControlPts',   UF.TGC,...
        'PulseInv',     UF.PulseInversion,...
        'TrigIn',       UF.TrigIn);
    %        'NbHostBuffer', 1);
    
catch ErrorMsg
    errordlg(ErrorMsg.message, ErrorMsg.identifier);
end

% ============================================================================ %
% ============================================================================ %

%% DO NOT CHANGE - Create and build the ultrasound sequence

try
    
    % Create the TPC object
    TPC = remote.tpc( ...
        'imgVoltage', ImagingVoltage, ...
        'imgCurrent', ImagingCurrent);
    
    % Create  and build the sequence
    Sequence           = usse.usse( 'TPC', TPC, 'acmo', Acmo, 'Popup',0 );
    [Sequence, NbAcqRx] = Sequence.buildRemote();
    
catch ErrorMsg
    errordlg(ErrorMsg.message, ErrorMsg.identifier);
end

%% beamforming

% ============================================================================ %
% BEAMFORMING
% ============================================================================ %

RxIds       = unique(Sequence.InfoStruct.mode(1).ModeRx(3,:));

bfStruct.fNumber = 1.2 ; % fnumber
bfStruct.linePitch = system.probe.Pitch*1e-3;
bfStruct.nbRecon = 192;
bfStruct.xOrigin = (system.probe.NbElemts - bfStruct.nbRecon)/2;
bfStruct.speedOfSound = common.constants.SoundSpeed;
bfStruct.piezoPitch = system.probe.Pitch*1e-3;
bfStruct.ProbeRadius = system.probe.Radius*1e-3;
bfStruct.demodFreq = UF.TxFreq*1e6;
bfStruct.rOrigin = (ImgInfo.Depth(1))*1e-3;
bfStruct.peakDelay = 0;
bfStruct.nbAngles = length(UF.xApex);
bfStruct.firstAcquiredChannel = Sequence.InfoStruct.rx(RxIds(1)).Channels(1)-1;
bfStruct.synthAcq = ( length(RxIds) == 2 );

bfStruct.lambda = bfStruct.speedOfSound/bfStruct.demodFreq;
bfStruct.pixelPitch = bfStruct.lambda/2;
bfStruct.nbPiezos = system.probe.NbElemts;
bfStruct.firstSample = Sequence.InfoStruct.event(end).skipSamples;
bfStruct.nbSamples = Sequence.InfoStruct.event(end).numSamples;
bfStruct.nbPixelsPerLine = round((ImgInfo.Depth(2)-ImgInfo.Depth(1))/bfStruct.pixelPitch*1e-3);
bfStruct.normMode=2;
bfStruct.frame_per_frame = -1;
bfStruct.idxTransmitToBeamform = 0;
bfStruct.usegpu = 1 ;
bfStruct.nbImages =UF.NbFrames;


for k=1:bfStruct.nbAngles
    bfStruct.xSources(k)   = UF.xApex(k)*1e-3;
    bfStruct.zSources(k)   = UF.zApex(k)*1e-3;
    bfStruct.timeOrigin(k) = double(Sequence.InfoStruct.tx(k).tof2Focus*1e-6);
end
bfStruct.channelOffset = Sequence.InfoStruct.mode(1).channelSize/2;

%%

clear parameters

parameters.NR = bfStruct.nbPixelsPerLine;
parameters.NTheta = bfStruct.nbRecon;

parameters.Nx = 512;
parameters.Nz = 512;

parameters.AngleOfFirstLine = (-system.probe.NbElemts / 2 +bfStruct.xOrigin)*bfStruct.piezoPitch/bfStruct.ProbeRadius;
parameters.AngleStep = bfStruct.linePitch/bfStruct.ProbeRadius;
parameters.FirstSampleAxialPosition = bfStruct.rOrigin + bfStruct.ProbeRadius;
parameters.RadialStep = bfStruct.pixelPitch;
parameters.OffSet = zeros(2,1, 'double');
parameters.SizeOfImg = zeros(2,1, 'double');
parameters.NbVois = 1;
parameters.SmoothDistance = 2;
parameters.SecondPass = 0;
parameters.OffSetX = 0;
parameters.OffSetZ = 0;
parameters.SizeOfImgX = 0;
parameters.SizeOfImgZ = 0;
parameters.ComputeOffSetAndSize = 1;
OUT = zeros(parameters.Nx,parameters.Nz,'single');

lutsize = zeros(1,2,'int32');

parameters_array = processing.SCV_Parser(parameters);
processing.ScanConversionTable(parameters_array,lutsize);

lutsize =lutsize(1);
%lutsize = int32(139924);

%%


parameters.OffSetX = parameters_array(10);
parameters.OffSetZ = parameters_array(11);
parameters.SizeOfImgX = parameters_array(12);
parameters.SizeOfImgZ = parameters_array(13);


NbVoistot = (parameters.NbVois*2)^2;

parameters.SecondPass = 1;
LUT_coeff = zeros(lutsize*NbVoistot,1, 'single');
LUT_X = zeros(lutsize,1, 'int16');
LUT_Z = zeros(lutsize,1, 'int16');
LUT_Vois1 = zeros(lutsize*NbVoistot,1, 'int16');
LUT_Vois2 = zeros(lutsize*NbVoistot,1, 'int16');

processing.ScanConversionTable(processing.SCV_Parser(parameters),LUT_coeff,LUT_X,LUT_Z,LUT_Vois1,LUT_Vois2,lutsize);


XAxis = parameters.OffSetX + (0:(parameters.Nx-1))*parameters.SizeOfImgX/(parameters.Nx-1);
ZAxis = parameters.OffSetZ + (0:(parameters.Nz-1))*parameters.SizeOfImgZ/(parameters.Nz-1);

%%


% figure();
% h=imagesc(XAxis,ZAxis,20*log10(abs(OUT/bfStruct.nbAngles)'));
% axis equal;
% axis tight;
% caxis([10 70]);
% colormap gray;
% drawnow;


%% Do NOT CHANGE - Sequence execution

Sequence = Sequence.initializeRemote('IPaddress', AixplorerIP,'InitSequence',0);

% Load sequence
Sequence = Sequence.loadSequence();

% Execute sequence
clear buff;
Sequence = Sequence.startSequence('Wait', 0);
zob = '0';

% Initialize remote on systems
buff=Sequence.getData('Realign', 0);

%% Beamforming

img_stack = {};
for k=1:repetition
    
    for nb = 1:bfStruct.nbAngles
        tic
        bfStruct.frame_per_frame = nb-1;
        IQ = beamformerSyntheticCurvedApex(bfStruct,buff.data,buff.alignedOffset);
        IQ2 = complex(IQ(1:2:end,:,:),IQ(2:2:end,:,:));
        img = abs(IQ2(:,:,1));
        
%         figure(2^(nb+1)-1); 
%         imagesc(20.*log10(double(img/bfStruct.nbAngles))); 
%         colormap gray;
%         caxis([20 80]);

        processing.ScanConvertLin(LUT_coeff,LUT_X,LUT_Z,LUT_Vois1,LUT_Vois2,parameters.NbVois,lutsize,single(img'),OUT);
        %set(h,'CData',20*log10(abs(OUT)'));
        figure(2^nb); cla;
        imagesc(XAxis,ZAxis,20*log10(abs(OUT/bfStruct.nbAngles)'));
        if nb>0
            hold on;
            plot(bfStruct.xSources(nb), bfStruct.zSources(nb), 'y.', ...
                'MarkerSize', 26);
            set(gca,'color','black')
        end
        axis equal;
        axis tight;
        caxis([20 80]);
        colormap gray;
        drawnow;
        zob = get(gcf,'currentcharacter') ;
        toc
        if(zob ~= '0')
            break ;
        end
    end
end

Sequence = Sequence.stopSequence('Wait', 0);

clearvars hdl;


%% Draw colorbox and characteristic points

for i = 1:1:length(UF.thApex);
    
    thetaS = UF.thApex(i);
    
    apx = tools.apex(box, probe, thetaS);
    
    % Drawings
    box.draw( 2^(i) );
    box.drawPoints( 2^(i) );
    apx.drawPoint( 2^(i) );
    probe.drawPoints( 2^(i) );
    apx.drawLines( box, probe, thetaS, 2^(i) );

end

                         
% %% Plot colorbox                                     
% 
% figure(2); hold on;
% if exist('hdl', 'var')
%     hdl(1) = box.draw(2, hdl(1));
% else
%     hdl(1) = box.draw(2);
% end
% 
% %% Plot delays
% 
% figure(42);
% nbApices = length(Sequence.InfoStruct.tx)
% cmap = hsv(nbApices)
% for i=1:nbApices
%     plot(Sequence.InfoStruct.tx(i, 1).Delays, 'o', 'Color', cmap(i, :));
%     hold on;
% end
% 

%% End of function

%end
