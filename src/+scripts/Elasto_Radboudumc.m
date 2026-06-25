% scripts.Elasto_Radboudumc
%   Build and run a sequence for the Ultra-COMPASS project
%
%   
%   with a push elusev and ultrafast imaging.
%
%   Prior running this script, please run 'addPathLegHAL V12'
%   When running this script, please select 'Add to Path' if asked
%
%   Note - This function is defined as a script of SCRIPTS package. It cannot be
%   used without the legHAL package developed by SuperSonic Imagine and without
%   a system with a REMOTE server running.
%
%   Copyright 2019 Supersonic Imagine

%% System parameters
AixplorerIP = '10.47.14.236'; % IP address of Aixplorer

%% Sequence parameters
Nb_loop = 3;    % Whole sequence repetition: integer [1 inf]
Trig_on_ECG = 1; % 0: start sequence immidiatly. 1: start sequence on ECG trig in signal

Nb_repeat_within_cardiac_cycle = 1;       % Sequence repetition after trig in signal [1 inf]
Pause_before_repeat = 10e-6; %250e-3;  % Pause after elasto part of the sequence [s]

ImagingVoltage = 60;           % imaging voltage [V]
ImagingCurrent = 2;            % security current limit [A]
PushVoltage    = 40;           % push voltage [V]
PushCurrent    = 10;           % security limit for push current [A]

% BMode plane wave ultrafast acquisition parameters
B.TxFreq       = 7.5;          % emission frequency [MHz]
B.NbHcycle     = 2;            % number of half cycles
B.FlatAngles   = sort( [0 -20 20 -15 15 10 -10 5 -5 2.5 -2.5] ); % flat angles [deg]
B.DutyCycle    = 1;            % duty cycle [0 1]
B.TXCenter     = system.probe.Pitch * system.probe.NbElemts/2;   % emission center [mm]
B.TXWidth      = system.probe.Pitch * system.probe.NbElemts;     % emission width [mm]
B.ApodFct      = 'none';       % emission apodization function (none, hanning)
B.RxCenter     = B.TXCenter ;  % reception center [mm]
B.RxWidth      = B.TXWidth;    % reception width [mm]
B.RxBandwidth  = 2;            % sampling mode (1 = 200%, 2 = 100%, 3 = 50%)
B.FIRBandwidth = 90;           % FIR receiving bandwidth [%] - center frequency = UFTxFreq
B.NbFrames     = 1;            % number of acquired images
B.PRF          = 8000;         % pulse frequency repetition [Hz] (0 for greatest possible)
B.TGC          = 100 * [ 0 3 5 7 8 8 8 8 ]; % TGC profile: vector of 8 values from surface to max depth, between 0 and 800


% ============================================================================ %
% Image parameters
ImgInfo.Depth(1)           = 1;   % initial depth [mm]
ImgInfo.Depth(2)           = 28;  % image depth [mm]
ImgInfo.CorrectFactorDepth = 0.4; % depth correction to estimate the acquisition duration

% ============================================================================ %
% Push parameters
PushInfo.TxFreq       = 5;             % push emission frequency [MHz]
PushInfo.DutyCycle    = 1;             % duty cycle [0, 1]
PushInfo.ApodFct      = 'none';        % apodization function
PushInfo.RatioFD      = 2;             % ratio F/D
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
UF.TxFreq       = 5;               % emission frequency [MHz]
UF.NbHcycle     = 6;               % # of half cycles
UF.FlatAngles   = [-4 -2 0 2 4];   % flat angles [deg]
UF.DutyCycle    = 1;               % duty cycle [0, 1]
UF.TXCenter     = system.probe.Pitch * system.probe.NbElemts/2; % emission center [mm]
UF.TXWidth      = system.probe.Pitch * system.probe.NbElemts;   % emission width [mm]
UF.ApodFct      = 'none';          % emission apodization function (none, hanning)
UF.RxCenter     = UF.TXCenter ;    % reception center [mm]
UF.RxWidth      = UF.TXWidth;      % reception width [mm]
UF.RxBandwidth  = 2;               % sampling mode (1 = 200%, 2 = 100%, 3 = 50%)
UF.FIRBandwidth = 60;              % FIR receiving bandwidth [%] - center frequency = UFTxFreq
UF.NbFrames     = 15;              % # of acquired images
UF.PRF          = 10000;           % pulse frequency repetition [Hz] (0 for greatest possible)
UF.TGC          = 800 * ones(1,8); % TGC profile


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
PRF = 2 * Nb_loop * Nb_repeat_within_cardiac_cycle * ...
    ( length(B.FlatAngles)  +  Nb_elasto_angles * UF.NbFrames ) ...
    / ( Nb_loop + (Nb_repeat_within_cardiac_cycle-1)*Pause_before_repeat )  ; % average PRF in Hz
    % Nb_loop considered to append every second (at BPM)
NbHalfCycle = B.NbHcycle; % HC

TxFreq = B.TxFreq ; % MHz
Impedance = 50 ; % ohm
DutyCycle = B.DutyCycle * 100 ;

Voltage = ImagingVoltage ; % V

alpha = 1.75 ;
K = 3.3e-7 ; % degrees/Joule

% NbFirings = length(B.FlatAngles) ;
% MaxSurfaceTemp = K * NbPiezo * NbFirings * (2*NbHalfCycle/TxFreq) * ...
%     (DutyCycle/100)^2 * Voltage^2 / Impedance  /  (NbFirings/PRF) * (TxFreq/TxFreq0)^alpha 

MaxSurfaceTemp_flats = K * NbPiezo * PRF * (2*NbHalfCycle/TxFreq) * ...
    (DutyCycle/100)^2 * Voltage^2 / Impedance * (TxFreq/TxFreq0)^alpha ;

%% push temperature elevation on probe estimation
TxFreq0 = 9 ; % MHz

NbPiezo = 50 ;
NbFirings = length(B.FlatAngles) ;
PRF = 2 * Nb_loop * Nb_repeat_within_cardiac_cycle * ...
    ( Nb_elasto_angles * length(PushInfo.PosZ) ) ...
    / ( Nb_loop + (Nb_repeat_within_cardiac_cycle-1)*Pause_before_repeat ) ; % average PRF in Hz
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
    error( [ 'Estimated temperature increase > 15 deg ' num2str(T_increase) ] )
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

    % Create the trig IN elusev
    E_trig_ECG = elusev.pause( ...
        'Pause', 5e-6, ...
        'TrigIn',Trig_on_ECG );

    % Create the pause elusev for sequence repetition
    E_pause_repeat = elusev.pause( ...
        'Pause', Pause_before_repeat );

    for n2 = 1:Nb_repeat_within_cardiac_cycle

        % Create the B mode acmo
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
            'ControlPts',   B.TGC, ...
            'NbHostBuffer', Nb_loop ); 

        % Create the SWE acmo: ultrafast acquisition and push elusev
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
            'elusev',       Push{n}, ...
            'NbHostBuffer', Nb_loop );
        end

        % Create the pause acmo
        AcmoList{end+1} = acmo.acmo( ...
            'elusev', E_pause_repeat );
    end

    AcmoList{1} = AcmoList{1}.setParam( 'elusev', E_trig_ECG );
        
catch ErrorMsg
    error(ErrorMsg.message, ErrorMsg.identifier);
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
    error(ErrorMsg.message, ErrorMsg.identifier);
end

% ============================================================================ %
% ============================================================================ %


%% Do NOT CHANGE - Sequence execution

%disp('press a key')
%pause

% Initialize remote on systems
Sequence = Sequence.initializeRemote('IPaddress', AixplorerIP, 'InitSequence',0 );

% Load sequence
Sequence = Sequence.loadSequence();

% Execute sequence
clear buffer;
Sequence = Sequence.startSequence('Wait', 1);

disp('Sequence started')

% wait end of sequence
pause(3*Nb_loop)

% Retrieve data
for k = 1:NbAcqRx * Nb_loop
    tic
    buffer{k} = Sequence.getData('Realign',2, 'Timeout',30); % Timeout in seconds
    data_time = toc;
    disp( [ 'Rate: ' num2str( length( buffer{k}.data )*2 / data_time / 1e6 ) ' MB/s' ] )
end

% Stop sequence
Sequence = Sequence.stopSequence('Wait', 1);

return

%% Plot push positions
figure(1)
hold off
for n=1:Nb_elasto_angles
    plot( PushInfo_angulated{n}.PosX, PushInfo_angulated{n}.PosZ, ...
        'xk', 'MarkerSize',12, 'LineWidth',2 )
    set(gca, 'Ydir','reverse')
    
    xlim( [ 0 system.probe.NbElemts*system.probe.Pitch ] )
    ylim( [ 0 25 ] )
    xlabel mm
    ylabel mm
    title( [ 'Push angle ' num2str(n) ] )
    pause
end

%% plot bufers

nbuf = 12;
nacq = size( buffer{nbuf}.RFdata, 3 )

figure(1)

for n = 1:nacq
    I = squeeze( buffer{nbuf}.RFdata(:,:,n) );
    imagesc( I )
    caxis( [-1 1]*1e3 )
    title( [ 'buffer ' num2str(nbuf) ', acq ' num2str(n) ] ) 
    pause(0.2)
end
