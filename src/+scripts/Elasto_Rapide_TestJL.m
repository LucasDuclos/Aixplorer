% scripts.Elasto
%   Build and run a sequence
%   with a push elusev and ultrafast imaging.
%
%   Prior running this script, please run 'addPathLegHAL V12'
%   When running this script, please select 'Add to Path' if asked
%
%   Note - This function is defined as a script of SCRIPTS package. It cannot be
%   used without the legHAL package developed by SuperSonic Imagine and without
%   a system with a REMOTE server running.

%% System parameters
AixplorerIP = '192.168.1.1'; % IP address of Aixplorer

%% Sequence parameters
Nb_loop     = 1;    % Whole sequence repetition: integer [1 inf] n epas dépasser 9
Trig_on_ECG = 0; % 0: start sequence immediatly. 1: start sequence on ECG trig in signal

Nb_repeat_within_cardiac_cycle = 5;       % Sequence repetition after trig in signal [1 inf] une sequence est constituée de 2 elasto donc 2 buffer ne pas dépasser 5
Pause_before_repeat            = 90e-3; %250e-3;  % Pause after elasto part of the sequence [s]

ImagingVoltage = 60;           % imaging voltage [V]
ImagingCurrent = 2;            % security current limit [A]
PushVoltage    = 40;           % push voltage [V]
PushCurrent    = 10;           % security limit for push current [A]

% BMode plane wave ultrafast acquisition parameters
B.TxFreq       = 7.5;          % emission frequency [MHz]
B.NbHcycle     = 2;            % number of half cycles
B.FlatAngles   = sort( [0 -1 1 -2 2 -3 3 -4 4 -5 5 -6 6 -7 7 -8 8 -9 9 -10 10 -11 11 -12 12 -13 13 -14 14 -15 15 -16 16 -17 17 -18 18 -19 19 -20 20] ); % flat angles [deg]
B.DutyCycle    = 1;            % duty cycle [0 1]
B.TXCenter     = 64*0.2; %system.probe.Pitch * system.probe.NbElemts/2;   % emission center [mm]
B.TXWidth      = 128*0.2; %system.probe.Pitch * system.probe.NbElemts;     % emission width [mm]
B.ApodFct      = 'none';       % emission apodization function (none, hanning)
B.RxCenter     = 64*0.2;  %B.TXCenter ;  % reception center [mm]
B.RxWidth      = 128*0.2; %B.TXWidth;    % reception width [mm]
B.RxBandwidth  = 2;            % sampling mode (1 = 200%, 2 = 100%, 3 = 50%)
B.FIRBandwidth = 90;           % FIR receiving bandwidth [%] - center frequency = UFTxFreq
B.NbFrames     = 1;            % number of acquired images
B.PRF          = 0;      %8000;         % pulse frequency repetition [Hz] (0 for greatest possible)
B.TGC          = 100 * [ 0 3 5 7 8 8 8 8 ]; % TGC profile: vector of 8 values from surface to max depth, between 0 and 800


% ============================================================================ %
% Image parameters
ImgInfo.Depth(1)           = 15;   % initial depth [mm]
ImgInfo.Depth(2)           = 30;  % image depth [mm]
ImgInfo.CorrectFactorDepth = 0.4; % depth correction to estimate the acquisition duration
ImgInfo.LinePitch          = 1 ;
ImgInfo.PitchPix           = 1 ;
ImgInfo.RxApod             = [0.4 0.6];
ImgInfo.NbColumnsPerPiezo  = 1;
ImgInfo.NbRowsPerLambda    = 1;

% Ultrafast image formation
UF.BeamformingKneeAngle = 0.30;
UF.BeamformingMaxAngle  = 0.40;
UF.DynRangedB           = 55;
UF.Threshold            = 100;

% ============================================================================ %
% Push parameters
PushInfo.TxFreq       = 5;             % push emission frequency [MHz]
PushInfo.DutyCycle    = 1;             % duty cycle [0, 1]
PushInfo.ApodFct      = 'none';        % apodization function    'hanning';
PushInfo.RatioFD      = 2;             % ratio F/D
PushInfo.PushDuration = 100;           % push duration [us]
PushInfo.Duration     = 125;           % event duration [us]
PushInfo.TxCenter     = 96*0.2; %system.probe.Pitch * system.probe.NbElemts/2;            % x-position of emission center [mm]
PushInfo.PosX         = 96*0.2; %system.probe.Pitch * system.probe.NbElemts/2;            % x-position of the focal point [mm]
PushInfo.PosZ         = [18 23 28];%30 35 40 45]; % z-position of the focal point [mm]


% Manage pushes

Nb_elasto_tirs      = 2;% 3 ;

% ============================================================================ %
% Ultrafast acquisition parameters
UF.TxFreq       = 7.5;               % emission frequency [MHz]
UF.NbHcycle     = 4;               % # of half cycles
UF.FlatAngles   = [-3 0 3];    %-4 -2 0 2 4];   % flat angles [deg]
UF.DutyCycle    = 1;               % duty cycle [0, 1]
UF.TXCenter     = 64*0.2;  %system.probe.Pitch * system.probe.NbElemts/2; % emission center [mm]
UF.TXWidth      = 128*0.2; %system.probe.Pitch * system.probe.NbElemts;   % emission width [mm]
UF.ApodFct      = 'none';          % emission apodization function (none, hanning)
UF.RxCenter     = 64*0.2; %UF.TXCenter ;    % reception center [mm]
UF.RxWidth      = 128*0.2;  %UF.TXWidth;      % reception width [mm]
UF.RxBandwidth  = 2;               % sampling mode (1 = 200%, 2 = 100%, 3 = 50%)
UF.FIRBandwidth = 60;              % FIR receiving bandwidth [%] - center frequency = UFTxFreq
UF.NbFrames     = 50;              % # of acquired images
UF.PRF          = 14000;           % pulse frequency repetition [Hz] (0 for greatest possible)
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
TxFreq0 = 9; % MHz

NbPiezo = 128 ;
PRF = 2 * Nb_loop * Nb_repeat_within_cardiac_cycle * ...
    ( length(B.FlatAngles)  +  Nb_elasto_tirs * UF.NbFrames ) ...
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
    ( Nb_elasto_tirs * length(PushInfo.PosZ) ) ...
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
    error( [ 'Estimated temperature increase > 15 °C ' num2str(T_increase) ] )
end

%% DO NOT CHANGE - Create acquisition mode
AcmoList = {};
Push = {};

try
    % Create the push elusev
    Push{1} = elusev.push( ...
        'TwFreq',       PushInfo.TxFreq, ...
        'DutyCycle',    PushInfo.DutyCycle, ...
        'ApodFct',      PushInfo.ApodFct, ...
        'RatioFD',      PushInfo.RatioFD, ...
        'PushDuration', PushInfo.PushDuration, ...
        'Duration',     PushInfo.Duration, ...
        'TxCenter',     PushInfo.TxCenter, ...
        'PosX',         PushInfo.PosX, ...
        'PosZ',         PushInfo.PosZ);
    
    Push{2} = elusev.push( ...
        'TwFreq',       PushInfo.TxFreq, ...
        'DutyCycle',    PushInfo.DutyCycle, ...
        'ApodFct',      PushInfo.ApodFct, ...
        'RatioFD',      PushInfo.RatioFD, ...
        'PushDuration', PushInfo.PushDuration, ...
        'Duration',     PushInfo.Duration, ...
        'TxCenter',     PushInfo.TxCenter - 64*0.2, ...
        'PosX',         PushInfo.PosX - 64*0.2, ...
        'PosZ',         PushInfo.PosZ);
    
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
        for n=1:Nb_elasto_tirs
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
            
            % Create the pause acmo
            AcmoList{end+1} = acmo.acmo( ...
                'elusev', E_pause_repeat );
            
        end
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
%pause(3*Nb_loop)

% Retrieve data
for k = 1:NbAcqRx * Nb_loop
    tic
    buffer{k} = Sequence.getData('Realign',2, 'Timeout',30); % Timeout in seconds
    data_time = toc;
    disp( [ 'Rate: ' num2str( length( buffer{k}.data )*2 / data_time / 1e6 ) ' MB/s' ] )
end

% Stop sequence
Sequence = Sequence.stopSequence('Wait', 1);

%return



%% beamforming
%% UF image formation

% Display parameters
XOrigin = double(min(Sequence.InfoStruct.rx(1).Channels - 1))* system.probe.Pitch; % reception left position
ZOrigin = ImgInfo.Depth(1);

% ============================================================================ %

% DO NOT CHANGE - Compute additional parameters
DeltaX = system.probe.Pitch / ImgInfo.NbColumnsPerPiezo;
DeltaZ = common.constants.SoundSpeed / 1000 /(UF.RxFreq / 4.0) / ImgInfo.NbRowsPerLambda;

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

for uu=1:length(buffer)
    tic
    [s,rf] = process_seq(Sequence,buffer{uu},ImgInfo,UF,XOrigin,ZOrigin);
    [IQ2(:,:,:,uu), velocity(:,:,:,uu)] = beamform_synth(s,rf,true,ImgInfo,PushInfo,1);
    disp(uu)
    toc
end


%% Plot push positions
% figure(1)
% hold off
% for n=1:Nb_elasto_tirs
%     plot( PushInfo_angulated{n}.PosX, PushInfo_angulated{n}.PosZ, ...
%         'xk', 'MarkerSize',12, 'LineWidth',2 )
%     set(gca, 'Ydir','reverse')
%
%     xlim( [ 0 system.probe.NbElemts*system.probe.Pitch ] )
%     ylim( [ 0 25 ] )
%     xlabel mm
%     ylabel mm
%     title( [ 'Push angle ' num2str(n) ] )
%     pause
% end

%% plot bufers

% nbuf = 6;
% nacq = size( buffer{nbuf}.RFdata, 3 )
%
% figure(1)
%
% for n = 1:nacq
%     I = squeeze( buffer{nbuf}.RFdata(:,:,n) );
%     imagesc( I )
%     caxis( [-1 1]*1e3 )
%     title( [ 'buffer ' num2str(nbuf) ', acq ' num2str(n) ] )
%     pause(0.2)
% end
