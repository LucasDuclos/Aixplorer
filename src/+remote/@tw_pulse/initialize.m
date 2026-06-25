function obj = initialize(obj, varargin)
   
% ============================================================================ %
% ============================================================================ %

current_class = 'remote.tw_pulse';

% Start error handling
try

% ============================================================================ %
% ============================================================================ %

% Initialize REMOTE.TW superclass
obj = initialize@remote.tw(obj, varargin{1:end});

% ============================================================================ %
% ============================================================================ %

%% Add new parameters

% Transmit frequency
Par = common.parameter( ...
    'TwFreq', ...
    'single', ...
    'sets the waveform frequency', ...
    {[0 system.hardware.MaxTxFreq]}, ...
    {'waveform frequency [MHz]'}, ...
    obj.Debug, current_class );
obj = obj.addParam(Par);

% ============================================================================ %

Par = common.parameter( ...
    'UseDualFreq', ...
    'int32', ...
    'add the 2*TwFreq emission frequency to transmitted waveform', ...
    {[0 1]}, ...
    {'0 or 1'}, ...
    obj.Debug, current_class );
Par = Par.setValue(0);
obj = obj.addParam(Par);

% ============================================================================ %

% # of half cycles (default = 2)
Par = common.parameter( ...
    'NbHcycle', ...
    'int32', ...
    'sets the number of half cycle', ...
    {[0 floor( double( system.hardware.MaxSamples ) / ( 0.5 * system.hardware.ClockFreq/2 /system.hardware.MaxTxFreq ) ) ]}, ...
    {'number of waveform half cycle'}, ...
    obj.Debug, current_class );
Par = Par.setValue(2);
obj = obj.addParam(Par);

% ============================================================================ %

% Waveform polarity (default = 1)
Par = common.parameter( ...
    'Polarity', ...
    'int32', ...
    'sets the waveform polarity', ...
    {[-1 1]}, ...
    {'-1: negative first arch, 1: positive first arch'}, ...
    obj.Debug, current_class );
Par = Par.setValue(1);
obj = obj.addParam(Par);

% ============================================================================ %

% Waveform sampling rate (default = 1)
Par = common.parameter( ...
    'txClock180MHz', ...
    'int32', ...
    'sets the waveform sampling rate', ...
    {0 1}, ...
    {'sampling rate = 90 MHz', 'sampling rate = 180 MHz'}, ...
    obj.Debug, current_class );
Par = Par.setValue(1);
obj = obj.addParam(Par);

% ============================================================================ %
% ============================================================================ %

%% End error handling
catch Exception
    
    % Exception in this method
    if ( isempty(Exception.identifier) )
        
        % Emit the new exception
        NewException = ...
            common.legHAL.GetException(Exception, class(obj), 'initialize');
        throw(NewException);

    % Re-emit previous exception
    else
        
        rethrow(Exception);
        
    end
    
end

% ============================================================================ %
% ============================================================================ %

end
