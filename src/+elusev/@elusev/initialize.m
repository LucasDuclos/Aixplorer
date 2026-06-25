function obj = initialize(obj, varargin)
   
% ============================================================================ %
% ============================================================================ %

current_class = 'elusev.elusev';

% Start error handling
try

% ============================================================================ %
% ============================================================================ %

% Initialize superclass
obj = initialize@common.remoteobj(obj, varargin{1:end});

% ============================================================================ %
% ============================================================================ %

%% Add new parameters

% Trigger in (default = 0)
Par = common.parameter( ...
    'TrigIn', ...
    'int32', ...
    'enables the trigger in', ...
    {0 1}, ...
    {'no trigger in', 'trigger in activated'}, ...
    obj.Debug, current_class );
Par = Par.setValue(0);
obj = obj.addParam(Par);

% ============================================================================ %

% Trigger out (default = 0)
Par = common.parameter( ...
    'TrigOut', ...
    'single', ...
    'enables the trigger out', ...
    {0 [1e-3 720.8]}, ...
    {'no trigger out', 'trigger duration [us]'}, ...
    obj.Debug, current_class );
Par = Par.setValue(0);
obj = obj.addParam(Par);

% ============================================================================ %

% Trigger out delay (default = 0)
Par = common.parameter( ...
    'TrigOutDelay', ...
    'single', ...
    'sets the trigger out delay', ...
    {[0 1000]}, ...
    {'trigger delay [us]'}, ...
    obj.Debug, current_class );
Par = Par.setValue(0);
obj = obj.addParam(Par);

% ============================================================================ %

% Triggers (in&out) on all events (default = 0)
Par = common.parameter( ...
    'TrigAll', ...
    'int32', ...
    'enables trriggers on all events', ...
    {0 1}, ...
    {'triggers on 1st event', 'triggers on all events'}, ...
    obj.Debug, current_class );
Par = Par.setValue(0);
obj = obj.addParam(Par);

% ============================================================================ %

% Number of repetition of the ELUSEV (default = 1)
Par = common.parameter( ...
    'Repeat', ...
    'int32', ...
    'sets the number of ELUSEV repetition', ...
    {[1 Inf]}, ...
    {'number of repetition'}, ...
    obj.Debug, current_class );
Par = Par.setValue(1);
obj = obj.addParam(Par);

% ============================================================================ %
% ============================================================================ %

%% Add new objects

% TX containers
obj = obj.addParam('tx', 'remote.tx');

% ============================================================================ %

% TW containers
obj = obj.addParam('tw', 'remote.tw');

% ============================================================================ %

% RX containers
obj = obj.addParam('rx', 'remote.rx');

% ============================================================================ %

% FC containers
obj = obj.addParam('fc', 'remote.fc');

% ============================================================================ %

% EVENT containers
obj = obj.addParam('event', 'remote.event');

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