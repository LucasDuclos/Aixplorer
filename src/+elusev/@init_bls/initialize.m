function obj = initialize(obj, varargin)
   
% ============================================================================ %
% ============================================================================ %

current_class = 'elusev.init_bls';

% Initialize superclass
try
    obj = initialize@elusev.elusev(obj, varargin{1:end});
    
    % Flat reception frequency
    Par = common.parameter( ...
        'RxFreq', ...
        'single', ...
        'sets the sampling frequency', ...
        {[1 60]}, ...
        {'flat reception frequency [MHz]'}, ...
        obj.Debug, current_class );
    obj = obj.addParam(Par);

catch exception
    uiwait(errordlg(exception.message, exception.identifier));
end

end