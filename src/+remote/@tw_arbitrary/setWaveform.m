% REMOTE.TW_ARBITRARY.SETWAVEFORM (PROTECTED)
%   Build the REMOTE.TW_ARBITRARY waveform.
%
%   WF = OBJ.SETWAVEFORM() returns the waveform WF of the REMOTE.TW_ARBITRARY
%   instance.
%
%   Note - This function is defined as a method of the remoteclass
%   REMOTE.TW_ARBITRARY. It cannot be used without all methods of the
%   remoteclass REMOTE.TW_ARBITRARY and all methods of its superclass REMOTE.TW
%   developed by SuperSonic Imagine and without a system with a REMOTE server
%   running.
%
%   Copyright 2010 Supersonic Imagine
%   Revision: 1.00 - Date: 2010/10/28

function varargout = setWaveform(obj, varargin)
   
% Start error handling
try

%% General controls on the method

% Check syntax
if nargout ~= 2
    
    % Build the prompt of the help dialog box
    ErrMsg = ['The ' upper(class(obj)) ' setWaveform function needs 2 ' ...
        'output argument.'];
    error(ErrMsg);
    
end

%% Build the waveform

% Retrieve parameters
TxElemts = obj.getParam('TxElemts');
Waveform = obj.getParam('Waveform');
RepeatCH = obj.getParam('RepeatCH');

% Control and build the final waveform
if RepeatCH
    Waveform_TxElemts = repmat(Waveform, [length(TxElemts), 1]);
else
    if length(TxElemts) == size(Waveform, 1)
        Waveform_TxElemts = Waveform;
    else
        % Build the prompt of the help dialog box
        ErrMsg = ['The arbitrary waveform is not defined for all transmit ' ...
            'elements'];
    end
end

Waveform_out = zeros(system.hardware.NbTxChan, size(Waveform,2), 'single');
Waveform_out( TxElemts, : ) = single(Waveform_TxElemts);

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
if ApodFct
    Waveform_out = obj.setApodization(Waveform_out, ApodFct, single(ones(1,length(TxElemts))), TxElemts);
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