function Delays = computeDelays(obj, FlatAngle)
% FlatAngle in rad

    % Initialize variables
    Delays     = zeros(1, system.hardware.NbTxChan);   %
    NbElemts   = system.probe.NbElemts;                %
    SoundSpeed = common.constants.SoundSpeed;          %

    % ============================================================================ %

    % Initialize variables
    PrPitch = system.probe.Pitch * 1e-3; % probe pitch [m]

    switch system.probe.Type

        case 'linear' % linear probe

            % Element cartesian coordinates
            ElemtXpos = ((1 : NbElemts) - 0.5) * PrPitch;     % element positions

            % Estimate the delays for FlatAngle
            Delays(1:NbElemts) = ElemtXpos * sin(FlatAngle) * 1e6 / SoundSpeed;

        case 'curved' % curved probe % control the flat law with Francois...

            % Element cartesian coordinates
            ElemtXpos = ((1 : NbElemts) - 0.5) * PrPitch;     % element positions

            % Estimate the delays for FlatAngle
            Delays(1:NbElemts) = ElemtXpos * sin(FlatAngle) * 1e6 / SoundSpeed;

        otherwise

            ErrMsg = ['The ' upper(class(obj)) ' setDelays function is not yet ' ...
                'implemented for ' lower(system.probe.Type) ' probes.'];
            error(ErrMsg);

    end

    % ============================================================================ %

    % Rescale delays
    Delays(1:NbElemts) = Delays(1:NbElemts) - min(Delays(1:NbElemts));

    % Set to 0 delays < TxClock
    Delays(Delays < 1/system.hardware.ClockFreq) = 0;
end