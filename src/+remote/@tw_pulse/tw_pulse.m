classdef tw_pulse < remote.tw
   
% ============================================================================ %
% ============================================================================ %

% Protected methods (accessible for subclasses)
methods ( Access = 'protected' )
    varargout = initalize(obj, varargin)   % build remoteclass REMOTE.TW_PULSE
    varargout = setWaveform(obj, varargin) % build the REMOTE.TW_PULSE waveform
    Waveform_Digital = dualFreqWfm( obj, Nsamples, NbHcycle, NDCycle)
end

% ============================================================================ %

% Class contructor
methods ( Access = 'public' )
    function obj = tw_pulse(varargin)
        
        % Label of the object
        if ( length(varargin) > 1 )
            if ( ~ischar(varargin{1}) || ~ischar(varargin{2}))
                varargin = { ...
                    'TW_PULSE', ...
                    'default periodic waveform', ...
                    varargin{1:end}};
            else
                varargin{1} = upper(varargin{1});
            end
        else
            varargin = { ...
                'TW_PULSE', ...
                'default periodic waveform', ...
                varargin{1:end}};
        end
        
        % Initialize object
        obj = obj@remote.tw(varargin{1:end});
        
    end
end

% ============================================================================ %
% ============================================================================ %

end