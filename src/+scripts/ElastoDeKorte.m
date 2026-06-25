% SCRIPTS.ELASTO (PUBLIC)
%   Build and run a sequence with a push elusev and ultrafast imaging.
%
%   Note - This function is defined as a script of SCRIPTS package. It cannot be
%   used without the legHAL package developed by SuperSonic Imagine and without
%   a system with a REMOTE server running.
%
%   Copyright 2010 Supersonic Imagine
%   Revision: 1.00 - Date: 2010/03/25

addpath /home/labo/git/rubi_v12/rubi/remoteClients/matlab/popeRubiBuild

%% Parameters definition

NThreads = 1 ;
Nb_loop = 3 ; % Whole sequence repetition

Nb_repeat_untrigged = 2 ;
Pause_before_repeat = 250e-3 ;

% System parameters
AixplorerIP = '192.168.3.58'; % IP address of Aixplorer
ImagingVoltage = 60;             % imaging voltage [V]
ImagingCurrent = 2;              % security current limit [A]
PushVoltage    = 40;             % push voltage [V]
PushCurrent    = 10;              % security limit for push current [A]

% BMode plane wave ultrafast acquisition parameters
B.TxFreq       = 7.5;            % emission frequency [MHz]
B.NbHcycle     = 2;               % # of half cycles
B.FlatAngles   = sort( [0 -20 20 -15 15 10 -10 5 -5 2.5 -2.5] ); % flat angles [deg]
B.DutyCycle    = 1;               % duty cycle [0, 1]
B.TXCenter     = system.probe.Pitch * system.probe.NbElemts/2;            % emission center [mm]
B.TXWidth      = system.probe.Pitch * system.probe.NbElemts;            % emission width [mm]
B.ApodFct      = 'none';          % emission apodization function (none, hanning)
B.RxCenter     = B.TXCenter ;            % reception center [mm]
B.RxWidth      = B.TXWidth;            % reception width [mm]
B.RxBandwidth  = 2;               % sampling mode (1 = 200%, 2 = 100%, 3 = 50%)
B.FIRBandwidth = 90;              % FIR receiving bandwidth [%] - center frequency = UFTxFreq
B.NbFrames     = 1;             % # of acquired images
B.PRF          = 8000;           % pulse frequency repetition [Hz] (0 for greatest possible)
B.TGC          = 100 * [ 0 3 5 7 8 8 8 8 ]; % TGC profile

% B image formation
B.BeamformingKneeAngle = 0.20;
B.BeamformingMaxAngle  = 0.30;
B.DynRangedB           = 55;
B.Threshold            = 100;

% ============================================================================ %
% Image parameters

ImgInfo.LinePitch = 1 ;
ImgInfo.PitchPix = 1 ;
ImgInfo.RxApod = [0.4 0.6];
ImgInfo.Depth(1)           = 1; % initial depth [mm]
ImgInfo.Depth(2)           = 28;  % image depth [mm]
ImgInfo.CorrectFactorDepth = 0.4; % depth correction to estimate the acquisition duration
ImgInfo.NbColumnsPerPiezo  = 4;
ImgInfo.NbRowsPerLambda    = 4;

% ============================================================================ %
% Push parameters

PushInfo.TxFreq       = 5 ;         % push emission frequency [MHz]
PushInfo.DutyCycle    = 1;           % duty cycle [0, 1]
PushInfo.ApodFct      = 'none';        % apodization function
PushInfo.RatioFD      = 2;           % ratio F/D
PushInfo.PushDuration = 150;           % push duration [us]
PushInfo.Duration     = 180;           % event duration [us]
PushInfo.PosZ         = (0:5:15) + 5 ; % z-position of the focal point [mm]

% Manage angulated pushes
Nb_elasto_angles = 3 ;
PushInfo_angulated{1} = PushInfo;
Push_X_on_probe       = system.probe.Pitch * 128 - 10 ;
Push_X_depest         = system.probe.Pitch * 192 - 10 ;
X_per_z = (Push_X_depest-Push_X_on_probe) / PushInfo_angulated{1}.PosZ(end) ;
for n = 1:length(PushInfo_angulated{1}.PosZ)
    z = PushInfo_angulated{1}.PosZ(n);
    PushInfo_angulated{1}.PosX(n)  = z * X_per_z + Push_X_on_probe ; % x-position of the focal point [mm]
end
PushInfo_angulated{1}.TxCenter     = Push_X_on_probe ;            % x-position of emission center [mm]

PushInfo_angulated{2} = PushInfo;
Push_X_on_probe       = system.probe.Pitch * 128 + 10 ;
Push_X_depest         = system.probe.Pitch * 64 + 10 ;
X_per_z = (Push_X_depest-Push_X_on_probe) / PushInfo_angulated{2}.PosZ(end) ;
for n = 1:length(PushInfo_angulated{2}.PosZ)
    z = PushInfo_angulated{2}.PosZ(n);
    PushInfo_angulated{2}.PosX(n)  = z * X_per_z + Push_X_on_probe ; % x-position of the focal point [mm]
end
PushInfo_angulated{2}.TxCenter     = Push_X_on_probe ;            % x-position of emission center [mm]

PushInfo_angulated{3} = PushInfo;
Push_X_on_probe       = system.probe.Pitch * 128 ;
Push_X_depest         = system.probe.Pitch * 128 ;
X_per_z = (Push_X_depest-Push_X_on_probe) / PushInfo_angulated{3}.PosZ(end) ;
for n = 1:length(PushInfo_angulated{3}.PosZ)
    z = PushInfo_angulated{3}.PosZ(n);
    PushInfo_angulated{3}.PosX(n)  = z * X_per_z + Push_X_on_probe ; % x-position of the focal point [mm]
end
PushInfo_angulated{3}.TxCenter     = Push_X_on_probe ;            % x-position of emission center [mm]

% ============================================================================ %

% Ultrafast acquisition parameters
UF.TxFreq       = 5;            % emission frequency [MHz]
UF.NbHcycle     = 6;               % # of half cycles
UF.FlatAngles   = [-4 -2 0 2 4];  % flat angles [deg]
UF.DutyCycle    = 1;               % duty cycle [0, 1]
UF.TXCenter     = system.probe.Pitch * system.probe.NbElemts/2;            % emission center [mm]
UF.TXWidth      = system.probe.Pitch * system.probe.NbElemts;            % emission width [mm]
UF.ApodFct      = 'none';          % emission apodization function (none, hanning)
UF.RxCenter     = UF.TXCenter ;            % reception center [mm]
UF.RxWidth      = UF.TXWidth;            % reception width [mm]
UF.RxBandwidth  = 2;               % sampling mode (1 = 200%, 2 = 100%, 3 = 50%)
UF.FIRBandwidth = 60;              % FIR receiving bandwidth [%] - center frequency = UFTxFreq
UF.NbFrames     = 15;             % # of acquired images
UF.PRF          = 10000;           % pulse frequency repetition [Hz] (0 for greatest possible)
UF.TGC          = 800 * ones(1,8); % TGC profile

% ============================================================================ %

% Ultrafast image formation
UF.BeamformingKneeAngle = 0.30;
UF.BeamformingMaxAngle  = 0.40;
UF.DynRangedB           = 55;
UF.Threshold            = 100;

% ============================================================================ %
% ============================================================================ %

%% DO NOT CHANGE - Additional parameters

% Additional parameters for ultrafast acquisition
UF.RxFreq = 4 * UF.TxFreq; % sampling frequency [MHz]

% Estimate RxDelay for ultrafast acquisition
UF.RxDelay = floor((1 - ImgInfo.CorrectFactorDepth) * ImgInfo.Depth(1) * 1e3 ...
    * UF.RxFreq / common.constants.SoundSpeed) / UF.RxFreq;

% Estimate RxDuration for ultrafast acquisition
UF.RxDuration = ceil((1 + ImgInfo.CorrectFactorDepth) ...
    * (2 * ImgInfo.Depth(2) - ImgInfo.Depth(1)) * 1e3 ...
    * (UF.RxFreq / 2^(UF.RxBandwidth - 1)) ...
    / (common.constants.SoundSpeed * 128)) * 128 ...
    / (UF.RxFreq / 2^(UF.RxBandwidth - 1));


% Additional parameters for ultrafast acquisition
B.RxFreq = 4 * B.TxFreq; % sampling frequency [MHz]

% Estimate RxDelay for ultrafast acquisition
B.RxDelay = floor((1 - ImgInfo.CorrectFactorDepth) * ImgInfo.Depth(1) * 1e3 ...
    * B.RxFreq / common.constants.SoundSpeed) / B.RxFreq;

% Estimate RxDuration for ultrafast acquisition
B.RxDuration = ceil((1 + ImgInfo.CorrectFactorDepth) ...
    * (2 * ImgInfo.Depth(2) - ImgInfo.Depth(1)) * 1e3 ...
    * (B.RxFreq / 2^(B.RxBandwidth - 1)) ...
    / (common.constants.SoundSpeed * 128)) * 128 ...
    / (B.RxFreq / 2^(B.RxBandwidth - 1));

% ============================================================================ %
% ============================================================================ %


%% flat temperature elevation on probe estimation
TxFreq0 = 9 ; % MHz

NbPiezo = 256 ;
NbFirings = length(B.FlatAngles) ;
PRF = 2 * Nb_loop * Nb_repeat_untrigged * ...
    ( length(B.FlatAngles)  +  Nb_elasto_angles * UF.NbFrames ) ...
    / ( Nb_loop-1 + (Nb_repeat_untrigged-1)*Pause_before_repeat )  ; % average PRF in Hz
    % Nb_loop considered to append every second (at BPM)
NbHalfCycle = B.NbHcycle; % HC

TxFreq = B.TxFreq ; % MHz
Impedance = 50 ; % ohm
DutyCycle = B.DutyCycle * 100 ;

Voltage = ImagingVoltage ; % V

alpha = 1.75 ;
K = 3.3e-7 ; % degrees/Joule

% MaxSurfaceTemp = K * NbPiezo * NbFirings * (2*NbHalfCycle/TxFreq) * ...
%     (DutyCycle/100)^2 * Voltage^2 / Impedance  /  (NbFirings/PRF) * (TxFreq/TxFreq0)^alpha 

MaxSurfaceTemp_flats = K * NbPiezo * PRF * (2*NbHalfCycle/TxFreq) * ...
    (DutyCycle/100)^2 * Voltage^2 / Impedance   * (TxFreq/TxFreq0)^alpha ;

%% push temperature elevation on probe estimation
TxFreq0 = 9 ; % MHz

NbPiezo = 50 ;
NbFirings = length(B.FlatAngles) ;
PRF = 2 * Nb_loop * Nb_repeat_untrigged * ...
    ( Nb_elasto_angles * length(PushInfo.PosZ) ) ...
    / ( Nb_loop-1 + (Nb_repeat_untrigged-1)*Pause_before_repeat ) ; % average PRF in Hz
    % Nb_loop considered to append every second (at BPM)
NbHalfCycle = PushInfo.Duration * PushInfo.TxFreq * 2 ; % HC

TxFreq = PushInfo.TxFreq ; % MHz
Impedance = 50 ; % ohm
DutyCycle = PushInfo.DutyCycle * 100 ;

Voltage = PushVoltage ; % V

alpha = 1.75 ;
K = 3.3e-7 ; % degrees/Joule

MaxSurfaceTemp_push = K * NbPiezo * PRF * (2*NbHalfCycle/TxFreq) * ...
    (DutyCycle/100)^2 * Voltage^2 / Impedance   * (TxFreq/TxFreq0)^alpha ;

%%
T_increase = MaxSurfaceTemp_flats + MaxSurfaceTemp_push ;

if T_increase > 15
    error( [ 'Estimated temperature increase > 15 deg' T_increase ] )
end

%% DO NOT CHANGE - Create acquisition mode
AcmoList = {};
Push = {};

try    
    % Create the push elusev
    for n=1:Nb_elasto_angles
        Push{n} = elusev.push( ...
        'TwFreq',       PushInfo_angulated{n}.TxFreq, ...
        'DutyCycle',    PushInfo_angulated{n}.DutyCycle, ...
        'ApodFct',      PushInfo_angulated{n}.ApodFct, ...
        'RatioFD',      PushInfo_angulated{n}.RatioFD, ...
        'PushDuration', PushInfo_angulated{n}.PushDuration, ...
        'Duration',     PushInfo_angulated{n}.Duration, ...
        'TxCenter',     PushInfo_angulated{n}.TxCenter, ...
        'PosX',         PushInfo_angulated{n}.PosX, ...
        'PosZ',         PushInfo_angulated{n}.PosZ);
    end

    E_trig_ECG = elusev.pause( ...
        'Pause', 5e-6, ...
        'TrigIn',1 );

    E_pause_repeat = elusev.pause( ...
        'Pause', Pause_before_repeat );

    for n2 = 1:Nb_repeat_untrigged

        AcmoList{end+1} = acmo.ultrafast( ...
            'TrigOut',      10, ...
            'TwFreq',       B.TxFreq, ...
            'NbHcycle',     B.NbHcycle, ...
            'FlatAngles',   B.FlatAngles, ...
            'DutyCycle',    B.DutyCycle, ...
            'TxCenter',     B.TXCenter, ...
            'TxWidth',      B.TXWidth, ...
            'ApodFct',      B.ApodFct, ...
            'RxFreq',       B.RxFreq, ...
            'RxDuration',   B.RxDuration, ...
            'RxDelay',      B.RxDelay, ...
            'RxCenter',     B.RxCenter, ...
            'RxWidth',      B.RxWidth, ...
            'RxBandwidth',  B.RxBandwidth, ...
            'FIRBandwidth', B.FIRBandwidth, ...
            'RepeatFlat',   B.NbFrames, ...
            'PRF',          B.PRF, ...
            'ControlPts',   B.TGC ); 

        % Create the ultrafast acquisition mode and add the push elusev
        for n=1:Nb_elasto_angles
            AcmoList{end+1} = acmo.ultrafast( ...
            'TrigOut',      10, ...
            'TwFreq',       UF.TxFreq, ...
            'NbHcycle',     UF.NbHcycle, ...
            'FlatAngles',   UF.FlatAngles, ...
            'DutyCycle',    UF.DutyCycle, ...
            'TxCenter',     UF.TXCenter, ...
            'TxWidth',      UF.TXWidth, ...
            'ApodFct',      UF.ApodFct, ...
            'RxFreq',       UF.RxFreq, ...
            'RxDuration',   UF.RxDuration, ...
            'RxDelay',      UF.RxDelay, ...
            'RxCenter',     UF.RxCenter, ...
            'RxWidth',      UF.RxWidth, ...
            'RxBandwidth',  UF.RxBandwidth, ...
            'FIRBandwidth', UF.FIRBandwidth, ...
            'RepeatFlat',   UF.NbFrames, ...
            'PRF',          UF.PRF, ...
            'ControlPts',   UF.TGC, ...
            'elusev',       Push{n} );
        end
        
        AcmoList{end+1} = acmo.acmo( ...
            'elusev', E_pause_repeat );
    end

    AcmoList{1} = AcmoList{1}.setParam( 'elusev', E_trig_ECG );
        
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
        'imgCurrent', ImagingCurrent, ...
        'pushVoltage', PushVoltage, ...
        'pushCurrent', PushCurrent);
    
    % Create  and build the sequence
    Sequence = usse.usse( 'TPC',TPC, ...
        'acmo', AcmoList, ...
        'Popup', 0, ...
        'Repeat', Nb_loop );
    [Sequence, NbAcqRx] = Sequence.buildRemote();
    
catch ErrorMsg
    errordlg(ErrorMsg.message, ErrorMsg.identifier);
end

% ============================================================================ %
% ============================================================================ %

%% Do NOT CHANGE - Sequence execution

% Initialize remote on systems
Sequence = Sequence.initializeRemote('IPaddress', AixplorerIP, 'InitSequence',0 );

% Load sequence
Sequence = Sequence.loadSequence();

% Execute sequence
clear buffer;
Sequence = Sequence.startSequence('Wait', 1);

% Retrieve data
for k = 1:NbAcqRx
    tic
    buffer{k} = Sequence.getData('Realign',2, 'Timeout',30); % Timeout in seconds
    data_time = toc;
    disp( [ 'Rate: ' num2str( length( buffer{k}.data )*2 / data_time / 1e6 ) ' MB/s' ] )
end
% Stop sequence
Sequence = Sequence.stopSequence('Wait', 1);

% return

% ============================================================================ %
% ============================================================================ %
%% B image

% DO NOT CHANGE - Image info
UFBfInfo.LinePitch = 1 / ImgInfo.NbColumnsPerPiezo;                       % dx / PiezoPitch
UFBfInfo.PitchPix  = 1 / ImgInfo.NbRowsPerLambda;                         % dz / lambda
UFBfInfo.Depth(1)  = ImgInfo.Depth(1);                                    % initial depth [mm]
UFBfInfo.Depth(2)  = ImgInfo.Depth(2);                                    % final depth [mm]
UFBfInfo.RxApod(1) = B.BeamformingKneeAngle * 4 * system.probe.Pitch ... % RxApodOne
    * 1000 / common.constants.SoundSpeed * B.RxFreq * 0.25;
UFBfInfo.RxApod(2) = B.BeamformingMaxAngle * 4 * system.probe.Pitch ...  % RxApodZero
    * 1000 / common.constants.SoundSpeed * B.RxFreq * 0.25;

UFBfInfo.RxApod(1) = min(2,UFBfInfo.RxApod(1));
UFBfInfo.RxApod(2) = min(2,UFBfInfo.RxApod(2));


% ============================================================================ %

% Beamform selected buffers
BFStruct_B = processing.lutSynthetic(Sequence, 1, UFBfInfo);
BFStruct_B.Info.DebugMode = 0 ;
BFStruct_B = processing.reconSynthFlat(BFStruct_B, buffer{1}, NThreads);

%% UF image formation

% Display parameters
disp( 'Beamforming B ...' )
XOrigin = double(min(Sequence.InfoStruct.rx(1).Channels - 1)) ... % reception left position
    * system.probe.Pitch;
ZOrigin = ImgInfo.Depth(1);

% ============================================================================ %

% DO NOT CHANGE - Compute additional parameters
DeltaX = system.probe.Pitch / ImgInfo.NbColumnsPerPiezo;
DeltaZ = common.constants.SoundSpeed / 1000 /(UF.RxFreq / 4.0) ...
    / ImgInfo.NbRowsPerLambda;

% DO NOT CHANGE - Image info
UFBfInfo.LinePitch = 1 / ImgInfo.NbColumnsPerPiezo;                       % dx / PiezoPitch
UFBfInfo.PitchPix  = 1 / ImgInfo.NbRowsPerLambda;                         % dz / lambda
UFBfInfo.Depth(1)  = ImgInfo.Depth(1);                                    % initial depth [mm]
UFBfInfo.Depth(2)  = ImgInfo.Depth(2);                                    % final depth [mm]
UFBfInfo.RxApod(1) = UF.BeamformingKneeAngle * 4 * system.probe.Pitch ... % RxApodOne
    * 1000 / common.constants.SoundSpeed * UF.RxFreq * 0.25;
UFBfInfo.RxApod(2) = UF.BeamformingMaxAngle * 4 * system.probe.Pitch ...  % RxApodZero
    * 1000 / common.constants.SoundSpeed * UF.RxFreq * 0.25;

UFBfInfo.RxApod(1) = min(2,UFBfInfo.RxApod(1));
UFBfInfo.RxApod(2) = min(2,UFBfInfo.RxApod(2));


% ============================================================================ %
BFStruct = processing.lutSynthetic(Sequence, 2, UFBfInfo);
BFStruct.Info.DebugMode = 0 ;

% Beamform selected buffers
disp( 'Beamforming Elasto ...' )
clear velocity
for n=1:Nb_elasto_angles
    BFStruct.ModeID = 1+n ;
    BFStruct = processing.reconSynthFlat(BFStruct, buffer{1+n}, NThreads);

    for k = 1:((BFStruct.Info.NImages)-1)
        correl = BFStruct.IQ(:,:,k).*conj(BFStruct.IQ(:,:,k+1));
        correl = filter2(ones(3*ImgInfo.NbRowsPerLambda,3*ImgInfo.NbColumnsPerPiezo)/9,correl);
        velocity{n}(:,:,k) = -angle(correl);
    end
end

% ============================================================================ %

% Compute X scales
UFNbBeamformedCols = size(BFStruct.IQ, 2);
UFXAxis            = XOrigin + DeltaX * (0 : (UFNbBeamformedCols - 1));

% Compute Z scales
UFNbBeamformedRows = size(BFStruct.IQ, 1);
UFZAxis            = ZOrigin + DeltaZ * (0 : (UFNbBeamformedRows - 1));


%% plot
IB = abs( BFStruct_B.IQ ) ;
figure( 1 )
imagesc( UFXAxis, UFZAxis, 20*log10( IB / max(IB(:)) ) )
caxis( [ -80 0] )
colormap gray
axis image
xlabel mm
ylabel mm

%%
for n=1:Nb_elasto_angles
    for k = 2:((BFStruct.Info.NImages)-1)
        figure( 2 )
        hold off
        imagesc( UFXAxis, UFZAxis, velocity{n}(:,:,k) )
        caxis( [ -pi pi ]/200 )
        colormap jet
        hold on
        plot( PushInfo_angulated{n}.PosX, PushInfo_angulated{n}.PosZ, ...
            'xk', 'MarkerSize',12, 'LineWidth',2 )
        axis image
        xlabel mm
        ylabel mm
        title( [ 'Push angle ' num2str(n) ] )
        pause(0.2)
    end
end
