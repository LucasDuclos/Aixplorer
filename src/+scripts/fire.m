function Wf1 = fire(voltage, nbHlfCycle, f0, dc, up, idx_pz_0, idx_pz_f)

% SCRIPTS.ARBITRARYACQ (PUBLIC)
%   Build and run a sequence with an ARBITRARY elusev.
%
%   Note - This function is defined as a script of SCRIPTS package. It cannot be
%   used without the legHAL package developed by SuperSonic Imagine and without
%   a system with a REMOTE server running.
%
%   Copyright 2010 Supersonic Imagine
%   Revision: 1.00 - Date: 2010/03/25

%% Parameters definition

%Arbitrary variables
% v1 = 20;
% f0      = 7.5;     % MHz
% dc = 0.5; up = 1;
% nbHlfCycle = 4;
% Pause   = 1 / f0; % us
% AcqDur  = 110; % us

% System parameters
AixplorerIP    = '192.168.3.131'; % IP address of the Aixplorer device '192.168.4.29'
ImagingVoltage = voltage;             % imaging voltage [V]
ImagingCurrent = 2;              % security current limit [A]

% ============================================================================ %

% Build the arbitrary waveform



% Create arbitrary waveform
SampFreq = system.hardware.ClockFreq;
NSamples = floor(SampFreq / (2*f0));
NUp = floor(dc*NSamples);

if up==1
    OneHalfCycle = [ones(1, NUp) zeros(1, NSamples-NUp)];
else
    OneHalfCycle = [-ones(1, NUp) zeros(1, NSamples-NUp)];
end

% Waveform
Wf1 = [];
for i = 1:1:nbHlfCycle
    Wf1 = horzcat(Wf1, (-1)^i*OneHalfCycle);
end

% Set waveform
Arbitrary.Waveform(:,:,1) = Wf1;
% Arbitrary.Waveform(:,:,2) = Wf2;
% Arbitrary.Waveform(:,:,3) = Wf3;

% ============================================================================ %

% DORT parameters
Arbitrary.Delays(:,1)  = 0+0*(idx_pz_0:idx_pz_f)'; % arbitrary delays [us]
% Arbitrary.Delays(:,2)  = 1+0*(1:system.probe.NbElemts)'; % arbitrary delays [us]
% Arbitrary.Delays(:,3)  = 2+0*(1:system.probe.NbElemts)'; % arbitrary delays [us]
Arbitrary.Pause        = 10;              % pause duration between events [us]
Arbitrary.PauseEnd     = 100;             % pause duration at the end of the elusev [us]
Arbitrary.ApodFct      = 'none';          % apodization function (none, bartlett, blackman,
                                          % connes, cosine, gaussian, hamming, hanning, welch)
Arbitrary.RxFreq       = 30;              % sampling frequency [MHz]
Arbitrary.RxDuration   = 98.2;            % acquisition duration [us]
Arbitrary.RxDelay      = 1;               % delay before acquisition [us]
Arbitrary.RxCenter     = 1;               % position of the receiving center [mm]
Arbitrary.RxWidth      = 1;               % width of the receiving window [mm]
Arbitrary.RxBandwidth  = 1;               % sampling mode (1 = 200%, 2 = 100%, 3 = 50%)
Arbitrary.FIRBandwidth = 90;              % FIR receiving bandwidth [%] - center frequency = UFTxFreq
Arbitrary.Repeat       = 1;               % repeat acquisitions
Arbitrary.TrigIn       = 0;               % enable trigger in
Arbitrary.TrigOut      = 1;               % duration of the trigger out [us]
Arbitrary.TrigOutDelay = 0;               % delay between the trigger out and the acquisition [us]
Arbitrary.TrigAll      = 0;               % enables triggers repeated for all events
Arbitrary.TGC          = 900 * ones(1,8); % TGC profile

% ============================================================================ %
% ============================================================================ %

%% DO NOT CHANGE - Create acquisition mode

try

    % Create the arbitrary elusev
    Elusev = elusev.arbitrary( ...
        'Waveform',     Arbitrary.Waveform, ...
        'Delays',       Arbitrary.Delays, ...
        'Pause',        Arbitrary.Pause, ...
        'PauseEnd',     Arbitrary.PauseEnd, ...
        'ApodFct',      Arbitrary.ApodFct, ...
        'RxFreq',       Arbitrary.RxFreq, ...
        'RxDuration',   Arbitrary.RxDuration, ...
        'RxDelay',      Arbitrary.RxDelay, ...
        'RxCenter',     Arbitrary.RxCenter, ...
        'RxWidth',      Arbitrary.RxWidth, ...
        'RxBandwidth',  Arbitrary.RxBandwidth, ...
        'FIRBandwidth', Arbitrary.FIRBandwidth, ...
        'TrigIn',       Arbitrary.TrigIn, ...
        'TrigOut',      Arbitrary.TrigOut, ...
        'TrigOutDelay', Arbitrary.TrigOutDelay, ...
        'TrigAll',      Arbitrary.TrigAll );

    % Create the ultrafast acquisition mode and add the arbitrary elusev
    Acmo = acmo.acmo( ...
        'ControlPts', Arbitrary.TGC, ...
        'elusev',     Elusev);

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
    Sequence           = usse.usse('TPC', TPC, 'acmo', Acmo, 'Popup',0 );
    [Sequence, ~] = Sequence.buildRemote();

catch ErrorMsg
    errordlg(ErrorMsg.message, ErrorMsg.identifier);
end

% return

% ============================================================================ %
% ============================================================================ %

%% Do NOT CHANGE - Sequence execution

try

    % Initialize remote on systems
    Sequence = Sequence.initializeRemote('IPaddress', AixplorerIP);

    % Load sequence
    Sequence = Sequence.loadSequence();

    % Execute sequence
    Sequence = Sequence.startSequence('Wait', 1);

    % Stop sequence
    Sequence = Sequence.stopSequence('Wait', 1);

catch ErrorMsg
    errordlg(ErrorMsg.message, ErrorMsg.identifier);
end

% ============================================================================ %
% ============================================================================ %


end

