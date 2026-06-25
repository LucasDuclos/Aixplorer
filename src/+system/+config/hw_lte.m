% SYSTEM.HARDWARE
%   Class containing hardware constants (Liver Therapy).
%
%   Constants:
%     - NAME corresponds to the name of the hardware.
%     - NBDAB corresponds to the number of DAB.
%     - NBADB corresponds to the number of ADB.
%     - NBTXPADB corresponds to the number of TX channels per ADBs.
%     - NBRXPADB corresponds to the number of RX channels per ADBs.
%     - NBTXCHAN corresponds to the number of TX channels.
%     - NBRXCHAN corresponds to the number of RX channels.
%     - CLOCKFREQ corresponds to the clock frequency.
%     - ADFILTER corresponds to the accepted decimation rates for the AD input
%       filter.
%     - ADRATE corresponds to the accepted AD sampling frequencies.
%     - MINNNOOP corresponds to the minimal duration of the noop variable.
%     - MAXNSAMPLES corresponds to the maximum number of samples per waveform.
%     - SGLsize corresponds to the maximum number of combinations between local
%       and host buffers.
%     - MAXEVENTS corresponds to the maximum number of events in a sequence.
%
%   Note - This class is defined as a member of the legHAL interface. It cannot
%   be used without a system with a REMOTE SERVER running.
%
%   Copyright 2010 Supersonic Imagine
%   Revision: 1.00 - Date: 2010/07/29

classdef hardware
   
% ============================================================================ %
% ============================================================================ %

% Hardware description
properties ( Constant = true )
    
    % Description
    Name = 'Liver Therapy';
    Tag  = 'LTE';
    
    % General constants
    NbDAB    = int32(2);   % number of DABs
    NbADB    = int32(16);  % number of ADBs
    NbTXpADB = int32(8);   % number of TX channels per ADB
    NbRXpADB = int32(4);   % number of RX channels per ADB
    NbTxChan = int32(256); % number of TX channels
    NbRxChan = int32(128); % number of RX channels

    % TX constants
    MaxTxFreq = single(4); % MHz

    % RX constants
    ClockFreq = single(180);    % clock frequency [MHz]
    ADFilter  = single(1:4);    % decimation rate for AD input filter
    ADRate    = single(3:9)'; % AD sampling frequency

    % Digital constants
    SGLsize  = int32(1024);  % number of combinations between local and host buffers
    MaxLocalBuffersSize = 1400e6; % in theory 2 GB avaliable on DAB but it crashed if higher than this

end
    
% ============================================================================ %

% Sequence constants
properties ( Constant = true )
    
    MinNoop    = single(5);     % minimal noop duration [us]
    MaxNoop    = single(5e5);   % maximal noop duration [us], in remote doc: 512*5200*0.2 = 532480 us
    MaxSamples = int32(4072);  % maximum number of samples for waveform
    MaxEvents  = int32(1e5); % maximum number of events in a sequence
    
end
    
% ============================================================================ %
% ============================================================================ %

end