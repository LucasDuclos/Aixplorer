%
% Generate WF as a successions of flats from angle min to max, with a
% settable minimal delay between two flats
% 
% angle is for all elements

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
TwFreq   = obj.getParam('TwFreq');
% UseDualFreq = obj.getParam('UseDualFreq');
DutyCycle = obj.getParam('DutyCycle');
NbHcycle = double(obj.getParam('NbHcycle'));
% SampFreq = 90 + 90 * double(obj.getParam('txClock180MHz'));
% Polarity = double(obj.getParam('Polarity'));
FlatAngles = obj.getParam('FlatAngles');
AngleId = obj.getParam('AngleId');

SampFreq = 180;
NbFlatAngles = length(FlatAngles) ;

% ============================================================================ %
% check parameters

if NbHcycle > floor( double( system.hardware.MaxSamples ) / ( 0.5 * system.hardware.ClockFreq/2 /TwFreq ) )
    error( [ 'NbHcycle is ' num2str(NbHcycle) ' but should be <= ' num2str( floor( double( system.hardware.MaxSamples ) / ( 0.5 * system.hardware.ClockFreq/2 /TwFreq ) ) ) ] )
end

if NbFlatAngles ~= 2^round(log2(NbFlatAngles))  &&  NbFlatAngles/12 ~= 2^round(log2(round(NbFlatAngles/12)))  &&  NbFlatAngles/20 ~= 2^round(log2(round(NbFlatAngles/20)))
    error( [ 'The number of angles N must satisfy that N, N/12 or N/20 is a power of 2, here N = ' num2str(NbFlatAngles) ] )
end

% ============================================================================ %

H = hadamard(NbFlatAngles);

% Build the temporal waveform without angulation
Nsamples = floor(0.5 * SampFreq / TwFreq);

% TODO : use DutyCycle
Waveform_Cos = cos(ceil((1 : Nsamples * NbHcycle) / Nsamples) * pi);
Waveform_Cos_length = length(Waveform_Cos);

%% compute delays
% Initialize variables
Delays     = zeros(1, system.hardware.NbTxChan);   %
NbElemts   = system.probe.NbElemts;                %
SoundSpeed = common.constants.SoundSpeed;          %
min_delay_between_flats = 4; % samples

% ============================================================================ %

% Initialize variables
PrPitch = system.probe.Pitch * 1e-3; % probe pitch [m]

FlatAngles_rad = FlatAngles * pi / 180; % flat angle [rad]

for n = 1:NbFlatAngles
    FlatAngle = FlatAngles_rad(n);

    switch system.probe.Type

        case 'linear' % linear probe

            % Element cartesian coordinates
            ElemtXpos = ((1 : NbElemts) - 0.5) * PrPitch;     % element positions

            % Estimate the delays for FlatAngle
            Delays(n, 1:NbElemts) = ElemtXpos * sin(FlatAngle) * 1e6 / SoundSpeed;

        case 'curved' % curved probe % control the flat law with Francois...

            % Element cartesian coordinates
            ElemtXpos = ((1 : NbElemts) - 0.5) * PrPitch;     % element positions

            % Estimate the delays for FlatAngle
            Delays(n, 1:NbElemts) = ElemtXpos * sin(FlatAngle) * 1e6 / SoundSpeed;

        otherwise

            ErrMsg = ['The ' upper(class(obj)) ' setDelays function is not yet ' ...
                'implemented for ' lower(system.probe.Type) ' probes.'];
            error(ErrMsg);
    end

    % Rescale delays
    Delays(n, 1:NbElemts) = Delays(n, 1:NbElemts) - min(Delays(n, 1:NbElemts));

    % Set to 0 delays < TxClock
    Delays(n, Delays(n,:) < 1/system.hardware.ClockFreq) = 0;
end

Delays_in_samples = round( Delays * SampFreq );

% minimise delta on channel NbElemts (assume angles are monotonaly increasing)
Delays_in_samples_optimal = 0 * Delays_in_samples ;
for n = 1:NbFlatAngles
    Delays_in_samples_optimal(n,:) = Delays_in_samples(n,:) - Delays_in_samples(n,NbElemts) ;
end

for n = 1:NbFlatAngles
    Delays_in_samples_absolute(n,:) = Delays_in_samples_optimal(n,:) + (NbFlatAngles-n) * (min_delay_between_flats + Waveform_Cos_length) ;
end
Delays_in_samples_absolute = Delays_in_samples_absolute - min(Delays_in_samples_absolute(:)) ;

obj.TWdelays = Delays_in_samples_absolute ; % used to beamform images

%% build waveform

% Build the initial full waveform
TotalNbSamples = max(Delays_in_samples_absolute(:)) + Waveform_Cos_length ;

if TotalNbSamples > system.hardware.MaxSamples
    error( [ 'Total number of samples is ' num2str(TotalNbSamples) ' but it must be <= system.hardware.MaxSamples which is ' num2str( system.hardware.MaxSamples ) ] )
end

if max( TxElemts ) > system.hardware.NbTxChan
    error( [ 'Emmission on element ' num2str(max( TxElemts )) ' but system.hardware.NbTxChan is ' num2str( system.hardware.NbTxChan ) ] )
end

% if length( Polarity ) == 1
%     Waveform = Polarity * Waveform_Cos;
%     Waveform_out = - repmat(Waveform, [length(TxElemts), 1]);
% 
% elseif length( Polarity ) == length( TxElemts )
%     if size(Polarity,1) == 1
%         Polarity = Polarity';
%     end
% 
%     Waveform = repmat(Waveform_Cos, [length(TxElemts), 1]);
%     Waveform_out = - repmat( Polarity, 1, length(Waveform_Cos) ) .* Waveform;
% 
% else
%     % Build the prompt of the help dialog box
%     ErrMsg = ['The length of the Polarity parameter for the class ' upper(class(obj)) ' setWaveform function ' ...
%         'must be 1 or equal to the number of transmitting elements.' ];
%     error(ErrMsg);
% 
% end

Waveform_out = single ( zeros(system.hardware.NbTxChan, TotalNbSamples) );

D = Delays_in_samples_absolute;
for n = 1:NbFlatAngles
    for e = TxElemts % 1:NbElemts
        Waveform_out( e, D(n,e)+1 : D(n,e)+Waveform_Cos_length ) = H(AngleId,n) * Waveform_Cos;
    end
end

% ============================================================================ %

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
