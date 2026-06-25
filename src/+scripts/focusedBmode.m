% clear all ; clear classes ;

% Sequence = usse.usse; Sequence = Sequence.selectProbe

% general parameters
NThreads = 12;
Debug = 0;

IPaddress = '192.168.3.140'; % IP address of Aixplorer
IPaddress     = 'local'; 

FrameRate  = 12; % Hz  (approx, SC6-1: 11 Hz because of scanconversion, L10_2: 15 Hz, SL15_4: 12 Hz)
repetition = FrameRate*10 ; % repeat sequence [1 inf]
DropFrames = 1;

B_repeat = 1 ; %
NbHostBuffer = 2;

% Voltage parameters
ImgVoltage  = 60;             % imaging voltage 0V]
ImgCurrent  = 0.1;            % security limit for imaging current [A]
PushVoltage = 10;             % imaging voltage [V]
PushCurrent = 5;              % security limit for imaging current [A]


% Image info
% fondamental
ImgInfo.TxPolarity     = [ 0 ]; % [ 0 1 ] : 0: Standard polarity, 1:Inverted polarity
ImgInfo.TxElemsPattern = [ 0 ]; % [ 0 2 ] : 0: All elements, 1: Odd elements, 2: Even elements
% puls inv (PI)
% ImgInfo.TxPolarity     = [ 0 1 ]; % [ 0 1 ] : 0: Standard polarity, 1:Inverted polarity
% ImgInfo.TxElemsPattern = [ 0 0 ]; % [ 0 2 ] : 0: All elements, 1: Odd elements, 2: Even elements
% checker board fondamental (PM)
% ImgInfo.TxPolarity     = [ 0 0 ]; % [ 0 1 ] : 0: Standard polarity, 1:Inverted polarity
% ImgInfo.TxElemsPattern = [ 1 2 ]; % [ 0 2 ] : 0: All elements, 1: Odd elements, 2: Even elements
% checker board + Pulse inversion (PMPÏ), only valid with RxBandwidth = 1
% ImgInfo.TxPolarity     = [ 0 1 0 ]; % [ 0 1 ] : 0: Standard polarity, 1:Inverted polarity
% ImgInfo.TxElemsPattern = [ 1 0 2 ]; % [ 0 2 ] : 0: All elements, 1: Odd elements, 2: Even elements
% ImgInfo.TxPolarity     = [ 0 0 1 1 ]; % [ 0 1 ] : 0: Standard polarity, 1:Inverted polarity
% ImgInfo.TxElemsPattern = [ 1 2 1 2 ]; % [ 0 2 ] : 0: All elements, 1: Odd elements, 2: Even elements

UseDualFreq = 0; % 1 to use diff THI, else 0

ImgInfo.SteeringAngle  = 0; % steering angle [not implemented yet]
ImgInfo.XOrigin        = 0 ; % xmin [mm]
ImgInfo.FirstLinePos   = 1 ; % position of the first line in elt [>=1]
ImgInfo.txPitch        = 1 ; % txpitch [elements]
ImgInfo.nbTx           = 192 ; % txpitch [elements]
ImgInfo.RxLinesPerTx   = 1;  % number of received line per tx (beamforming)
ImgInfo.PitchPix       = 0.5;  % dz / lambda
ImgInfo.Depth(1)       = 1;   % initial depth [mm]
ImgInfo.Depth(2)       = 80;  % final depth [mm]
ImgInfo.RxApod(1)      = 0.4; % RxApodOne
ImgInfo.RxApod(2)      = 0.6; % RxApodZero

ImgInfo.pulseInversion = length( unique( ImgInfo.TxPolarity ) ) - 1 ;

% Emission parameters
TwFreq          = 4.5/(1+ImgInfo.pulseInversion);      % transmit frequency [MHz]
NbHcycle        = 4;         % nb tx half cycle
ApodFct         = 'hanning';    % Apod function on the transmit aperture
TxFDRatio       = 3.5;         % F/D ratio (transmit)
DutyCycle       = 1;         % tx Duty cycle [0-1]
Focus           = prof/2 ;       % focal depth [mm]

% ============================================================================ %

% Receive parameters
RxBandwidth  = 2;          % 1 = 200%, 2 = 100%
RxFilterBandwidth = 90 ;
TGC = [1 400 600 600 600 600];

%% converting ImgInfo parameters into acmo parameters

RxFreq       = 4 * TwFreq * (1+ImgInfo.pulseInversion - .5 * UseDualFreq); % MHz

NumSamples   = ceil(1.4 * (2 * ImgInfo.Depth(2) - ImgInfo.Depth(1)) * 1e3 ...
    * RxFreq / (common.constants.SoundSpeed * 2^(RxBandwidth-1) * 128)) * 128;
SkipSamples  = floor(0. * ImgInfo.Depth(1) * RxFreq * 1e3 ...
    / common.constants.SoundSpeed);

RxDuration = NumSamples / RxFreq * 2^(RxBandwidth-1);
RxDelay    = SkipSamples / RxFreq;

%% built sequence

Acmo = acmo.bfoc( ...
    'TwFreq',       TwFreq, ...
    'UseDualFreq',  UseDualFreq, ...
    'NbHcycle',     NbHcycle, ...
    'SteerAngle',   ImgInfo.SteeringAngle, ...
    'Focus',        Focus, ...
    'PRF',          FrameRate, ...
    'DutyCycle',    DutyCycle, ...
    'XOrigin',      ImgInfo.XOrigin, ...
    'FirstLinePos', ImgInfo.FirstLinePos,...
    'NbFocLines',   ImgInfo.nbTx,...
    'TxLinePitch',  ImgInfo.txPitch,...
    'FDRatio',      TxFDRatio, ...
    'ApodFct',      ApodFct, ...
    'RxFreq',       RxFreq, ...
    'RxDuration',   RxDuration, ...
    'RxDelay',      RxDelay, ...
    'RxBandwidth',  RxBandwidth, ...
    'FIRBandwidth', RxFilterBandwidth, ...
    'PulseInv',     ImgInfo.pulseInversion, ...
    'TxPolarity',   ImgInfo.TxPolarity, ...
    'TxElemsPattern', ImgInfo.TxElemsPattern, ...
    'ControlPts',   TGC, ...
    'Repeat', B_repeat, ...
    'NbHostBuffer', NbHostBuffer, ...
    0);

TPC = remote.tpc( ...
    'imgVoltage',   ImgVoltage, ...
    'imgCurrent',   ImgCurrent, ...
    'pushVoltage',  PushVoltage, ...
    'pushCurrent',  PushCurrent, ...
    0);

% add a bfoc acmo and built the sequence
bfocId = 1;
Sequence = usse.usse( 'acmo',Acmo, 'TPC',TPC, 'Popup',0 );
Sequence = Sequence.setParam( 'Repeat',repetition, 'DropFrames',DropFrames );

[Sequence, NbAcq] = Sequence.buildRemote;

% Sequence.plot_diagram; xlim( [ 0 20e5 ] ); return

%% build beamforming luts if needed
BFStruct.Info.DebugMode = 0 ;
disp('LUT computation');
BFStruct = processing.lutBfoc( Sequence, bfocId, ImgInfo, NThreads );
if system.probe.Radius > 0
    Data = processing.lutScanConvertImage( BFStruct, ImgInfo, 1024, 1024 ); % TODO: remove hardcoded size
end

% clear buffer if necessary
clear buffer

%% Setup RF buffer
all_bufs = {};
buf_size = Sequence.InfoStruct.mode(1).channelSize(1) * system.hardware.NbRxChan / 2 ;
all_bufs{1} = int16( zeros( 1, buf_size, 'int16' ) );

BFStruct.Info.DebugMode = 0 ;

%   Initialize & load sequence remote
Sequence = Sequence.initializeRemote( 'IPaddress',IPaddress, 'InitSequence',1 );
% rmpath /home/labo/git_master/us/legHAL/src/../lib/libRemoteV9
% addpath ../../rubi/remoteClients/matlab/remoteBuild/

Sequence = Sequence.loadSequence();

% Start sequence
disp('start');
Sequence = Sequence.startSequence();

clear buf
time_to_get_data = zeros( 1, NbAcq*repetition );
time_to_loop = time_to_get_data;

for k=1:NbAcq*repetition
    k
    tic_loop = tic;
    disp('Getting buffer');

    buffer_data = all_bufs{1};

    tic_get_data = tic;
    buf{k} = Sequence.getData('Realign',0, 'Timeout',120, 'Buffer',buffer_data, 'DataFormat',Sequence.getParam('DataFormat'), 'Debug',Debug );
    time_to_get_data(k) = toc(tic_get_data);

    % Stop sequence
    disp('reconstruction');

    tic_recon = tic;
    BFStruct = processing.reconBFoc( BFStruct, buf{k} );
    time_to_recon = toc(tic_recon);

    tic_scan_conv = tic;
    if system.probe.Radius > 0
        Data = processing.ScanConvertImage( Data, abs(BFStruct.IQ) );
        img = Data.img ;
        X = Data.XAxis;
        Z = Data.ZAxis;
    else
        X = ImgInfo.XOrigin  + ( 0:(ImgInfo.nbTx * ImgInfo.RxLinesPerTx - 1) ) * ( system.probe.Pitch * ImgInfo.txPitch / ImgInfo.RxLinesPerTx ) ;
        Z = ImgInfo.Depth(1) + ( 0:(BFStruct.Lut(1).NPixLine - 1 ) ) * 1e-3 * common.constants.SoundSpeed / (RxFreq/4) * ImgInfo.PitchPix ;
        img = 20*log10( abs( BFStruct.IQ ) );
    end
    time_to_scan_conv = toc(tic_scan_conv);

    tic_plot = tic;
    if k == 1
        figure(1);
        h = imagesc(X,Z,img(:,:,1));

        axis equal ;
        axis tight ;
        caxis ([70 140]);

        colormap gray ;
    else
        set(h,'CData', img(:,:,1) )
    end
    drawnow
    time_to_plot = toc(tic_plot);

    time_to_loop(k) = toc(tic_loop);
    
    if Debug
        disp( [ 'Total time to loop : ' num2str(time_to_loop(k)) ] )
    end
end
Sequence = Sequence.stopSequence('Wait', 0);

if Debug
    fr_recon = 1 / time_to_recon

    BW_MB_s = length(buf{1}.data)*2 ./ time_to_get_data ./1024./1024 % MB / s
    fr_get_data = 1 ./ time_to_get_data
    fr_loop = 1 ./ time_to_loop
    fr_loop_mean = mean(fr_loop)
end

return

if ImgInfo.pulseInversion
    ch = size( buffer.data, 2) / 2;
    win_sample = 1:1000;
    acq_offset = 128;

    display_range = [-1 1]*800;

    figure; imagesc( squeeze( buffer.data(:,:,acq_offset) ) ); caxis( display_range )

    A = squeeze( buffer.data(win_sample,ch,acq_offset+1) ); 
    B = squeeze( buffer.data(win_sample,ch,acq_offset+2) );
    C = A+B;

    figure
    plot( [ A B C ] ); ylim( display_range )

    figure; spectrogram(double(A), 128,[],[], RxFreq ); caxis auto; the_caxis = caxis;
    figure; spectrogram(double(C), 128,[],[], RxFreq ); caxis( the_caxis )
end

