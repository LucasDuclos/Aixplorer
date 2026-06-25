% REMOTE.TW_PULSE.SETWAVEFORM (PROTECTED)
%   Build the REMOTE.TW_PULSE waveform.
%
%   WF = OBJ.SETWAVEFORM() returns the waveform WF of the REMOTE.TW_PULSE
%   instance.
%
%   Note - This function is defined as a method of the remoteclass
%   REMOTE.TW_PULSE. It cannot be used without all methods of the remoteclass
%   REMOTE.TW_PULSE and all methods of its superclass REMOTE.TW developed by
%   SuperSonic Imagine and without a system with a REMOTE server running.
%
%   Copyright 2010 Supersonic Imagine
%   Revision: 1.00 - Date: 2010/08/02

function varargout = setWaveform(obj, varargin)

% Start error handling
try

%% General controls on the method

% Check syntax
if ( nargout ~= 2 )

    % Build the prompt of the help dialog box
    ErrMsg = ['The ' upper(class(obj)) ' setWaveform function needs 2 ' ...
        'output argument.'];
    error(ErrMsg);

end

%% Build the waveform

% Retrieve parameters
TxElemts = obj.getParam('TxElemts');
TwFreq   = obj.getParam('TwFreq');
UseDualFreq = obj.getParam('UseDualFreq');
DutyCycle = obj.getParam('DutyCycle');
NbHcycle = double(obj.getParam('NbHcycle'));
SampFreq = 90 + 90 * double(obj.getParam('txClock180MHz'));
Polarity = double(obj.getParam('Polarity'));

% ============================================================================ %

if NbHcycle > floor( double( system.hardware.MaxSamples ) / ( 0.5 * system.hardware.ClockFreq/2 /TwFreq ) )
    error( [ 'NbHcycle is ' num2str(NbHcycle) ' but should be <= ' num2str( floor( double( system.hardware.MaxSamples ) / ( 0.5 * system.hardware.ClockFreq/2 /TwFreq ) ) ) ] )
end

% Build the temporal waveform
Nsamples     = floor(0.5 * SampFreq / TwFreq);

if UseDualFreq
    Waveform_Cos = obj.dualFreqWfm(Nsamples, NbHcycle, DutyCycle*Nsamples);
else
    Waveform_Cos = cos(ceil((1 : Nsamples * NbHcycle) / Nsamples) * pi);
end

if max( TxElemts ) > system.hardware.NbTxChan
    error( [ 'Emmission on element ' num2str(max( TxElemts )) ' but system.hardware.NbTxChan is ' num2str( system.hardware.NbTxChan ) ] )
end

if length( Polarity ) == 1
    Waveform = Polarity * Waveform_Cos;
    Waveform_TxElemts = - repmat(Waveform, [length(TxElemts), 1]);

elseif length( Polarity ) == length( TxElemts )
    if size(Polarity,1) == 1
        Polarity = Polarity';
    end

    Waveform = repmat(Waveform_Cos, [length(TxElemts), 1]);
    Waveform_TxElemts = - repmat( Polarity, 1, length(Waveform_Cos) ) .* Waveform;

else
    % Build the prompt of the help dialog box
    ErrMsg = ['The length of the Polarity parameter for the class ' upper(class(obj)) ' setWaveform function ' ...
        'must be 1 or equal to the number of transmitting elements.' ];
    error(ErrMsg);

end


% Update apodization window with DUTYCYCLE
if length(DutyCycle) == 1
    DutyCycle = DutyCycle .* ones(size(TxElemts));
elseif length(DutyCycle) ~= length(TxElemts)
    % Build the prompt of the help dialog box
    ErrMsg = [ 'The DUTYCYCLE is defined for ' num2str(length(DutyCycle)) ...
        'while there are ' num2str(length(TxElemts)) ' firing elements.' ];
    error(ErrMsg);
end

%min_delay_between_pulses = 4; % samples
Waveform_length = length(Waveform_Cos) ; %+ min_delay_between_pulses + size( Additional_pulse, 2 );

Waveform_out = zeros( system.hardware.NbTxChan, Waveform_length );
Waveform_out( TxElemts, 1:length(Waveform_Cos) ) = single(Waveform_TxElemts);

switch obj.getParam('ApodFct')
    case 'bartlett'
        ApodFct = 1;
    case 'blackman'
        ApodFct = 2;
    case 'connes'
        ApodFct = 3;
    case 'cosine'
        ApodFct = 4;
    case 'gaussian'
        ApodFct = 5;
    case 'hamming'
        ApodFct = 6;
    case 'hanning'
        ApodFct = 7;
    case 'welch'
        ApodFct = 8;
    otherwise
        ApodFct = 0;
end

% ApodFct = 6 ;
% DutyCycle = 0.9 ;

% best is -1:1, -5 dB @ 123 MHz, but peigne
% for k = -1:1
%     Waveform_out( :, 90+k:90:end ) = 0 ;
% end

if UseDualFreq == 0 % && ApodFct
    Waveform_out = obj.setApodization( single(Waveform_out), ApodFct, DutyCycle, TxElemts);
end

varargout{1} = obj;
varargout{2} = single( Waveform_out );


%% End error handling
catch Exception

    % Exception in this method
    if ( isempty(Exception.identifier) )

        % Emit the new exception
        NewException = ...
            common.legHAL.GetException(Exception, class(obj), 'setWaveform');
        throw(NewException);

    % Re-emit previous exception
    else

        rethrow(Exception);

    end

end

end
