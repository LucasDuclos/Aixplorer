function [ Z, noiseFloor, noiseStDev, brightness, noise_rf, brightness_rf ] = bmode_snr(ip, voltage, duty_cycle, demodFreq, filterBw, depth, vec_TGC, p, bfg)

    addpath /home/labo/git/us/legHAL/lib/popeRubi

    % general parameters
    NThreads = 12;
    Debug = 0;

    IPaddress = ip; % IP address of Aixplorer
    %IPaddress     = 'local'; 

    FrameRate  = 12; % Hz  (approx, SC6-1: 11 Hz because of scanconversion, L10_2: 15 Hz, SL15_4: 12 Hz)
    repetition = 1;%FrameRate*10 ; % repeat sequence [1 inf]
    DropFrames = 1;

    B_repeat = 1 ; %
    NbHostBuffer = 2;

    % Voltage parameters
    ImgVoltage  = voltage;             % imaging voltage 0V]
    ImgCurrent  = 0.1;            % security limit for imaging current [A]
    PushVoltage = 10;             % imaging voltage [V]
    PushCurrent = 5;              % security limit for imaging current [A]


    % Image info
    % fondamental
    ImgInfo.TxPolarity     = [ 0 ]; % [ 0 1 ] : 0: Standard polarity, 1:Inverted polarity
    ImgInfo.TxElemsPattern = [ 0 ]; % [ 0 2 ] : 0: All elements, 1: Odd elements, 2: Even elements
    % puls inv (PI)
    % ImgInfo.TxPolarity     = [ 0 1 ]; % [ 0 1 ] : 0: Standard polarity, 1:Inverted polarity
    % ImgInfo.TxElemsPattern = [ 0 0 ]; % [ 0 2 ] : 0: All elements, 1: Odd elements, 2: Even elements
    % checker board fondamental (PM)
    % ImgInfo.TxPolarity     = [ 0 0 ]; % [ 0 1 ] : 0: Standard polarity, 1:Inverted polarity
    % ImgInfo.TxElemsPattern = [ 1 2 ]; % [ 0 2 ] : 0: All elements, 1: Odd elements, 2: Even elements
    % checker board + Pulse inversion (PMPÏ), only valid with RxBandwidth = 1
    % ImgInfo.TxPolarity     = [ 0 1 0 ]; % [ 0 1 ] : 0: Standard polarity, 1:Inverted polarity
    % ImgInfo.TxElemsPattern = [ 1 0 2 ]; % [ 0 2 ] : 0: All elements, 1: Odd elements, 2: Even elements
    % ImgInfo.TxPolarity     = [ 0 0 1 1 ]; % [ 0 1 ] : 0: Standard polarity, 1:Inverted polarity
    % ImgInfo.TxElemsPattern = [ 1 2 1 2 ]; % [ 0 2 ] : 0: All elements, 1: Odd elements, 2: Even elements

    UseDualFreq = 0; % 1 to use diff THI, else 0
    steeringAngleList = p.steeringAngle;
   
    ImgInfo.SteeringAngle  = steeringAngleList; % steering angle [not implemented yet]
    ImgInfo.XOrigin        = 0 ; % xmin [mm]
    ImgInfo.FirstLinePos   = 1 ; % position of the first line in elt [>=1]
    ImgInfo.txPitch        = p.txLinePitch ; % txpitch [elements]
    ImgInfo.nbTx           = p.nbTXChannels ; % txpitch [elements]
    ImgInfo.RxLinesPerTx   = bfg.nbLines;  % number of received line per tx (beamforming)
    ImgInfo.PitchPix       = bfg.pixelPitch;  % dz / lambda
    ImgInfo.Depth(1)       = p.zOrigin;   % initial depth [mm]
    ImgInfo.Depth(2)       = depth;  % final depth [mm]
    ImgInfo.RxApod(2)      = 1.2; % RxApodZero < 2
    ImgInfo.RxApod(1)      = ImgInfo.RxApod(2)/8; % RxApodOne < 2
    

    ImgInfo.pulseInversion = length( unique( ImgInfo.TxPolarity ) ) - 1 ;

    % Emission parameters
    TwFreq          = p.txFreq;      % transmit frequency [MHz]
    NbHcycle        = p.nb_halfCycles;         % nb tx half cycle
    ApodFct         = 'hanning';    % Apod function on the transmit aperture
    TxFDRatio       = 4;         % F/D ratio (transmit)
    DutyCycle       = duty_cycle;         % tx Duty cycle [0-1]
    Focus           = p.focusMm ;       % focal depth [mm]

    % ============================================================================ %

    % Receive parameters
    RxBandwidth  = 2;          % 1 = 200%, 2 = 100%
    RxFilterBandwidth = filterBw ;
    TGC = vec_TGC;
   
    %% converting ImgInfo parameters into acmo parameters

    % RxFreq       = 4 * TwFreq * (1+ImgInfo.pulseInversion - .5 * UseDualFreq); % MHz
    RxFreq = 4*demodFreq;

    NumSamples   = ceil(1.4 * (2 * ImgInfo.Depth(2) - ImgInfo.Depth(1)) * 1e3 ...
        * RxFreq / (common.constants.SoundSpeed * 2^(RxBandwidth-1) * 128)) * 128;
    SkipSamples  = floor(0. * ImgInfo.Depth(1) * RxFreq * 1e3 ...
        / common.constants.SoundSpeed);

    RxDuration = NumSamples / RxFreq * 2^(RxBandwidth-1);
    RxDelay    = SkipSamples / RxFreq;

    for DutyCycleLoop = 0:1
        %% built sequence
        if DutyCycleLoop == 0
            DutyCycle = 0;
        else
            DutyCycle = duty_cycle;
        end
        Acmo = acmo.bfoc( ...
            'TwFreq',       TwFreq, ...
            'UseDualFreq',  UseDualFreq, ...
            'NbHcycle',     NbHcycle, ...
            'SteerAngle',   ImgInfo.SteeringAngle, ...
            'Focus',        Focus, ...
            'PRF',          FrameRate, ...
            'DutyCycle',    DutyCycle, ...
            'XOrigin',      ImgInfo.XOrigin, ...
            'FirstLinePos', ImgInfo.FirstLinePos,...
            'NbFocLines',   ImgInfo.nbTx,...
            'TxLinePitch',  ImgInfo.txPitch,...
            'FDRatio',      TxFDRatio, ...
            'ApodFct',      ApodFct, ...
            'RxFreq',       RxFreq, ...
            'RxDuration',   RxDuration, ...
            'RxDelay',      RxDelay, ...
            'RxBandwidth',  RxBandwidth, ...
            'FIRBandwidth', RxFilterBandwidth, ...
            'PulseInv',     ImgInfo.pulseInversion, ...
            'TxPolarity',   ImgInfo.TxPolarity, ...
            'TxElemsPattern', ImgInfo.TxElemsPattern, ...
            'ControlPts',   TGC, ...
            'Repeat', B_repeat, ...
            'NbHostBuffer', NbHostBuffer, ...
            0);

        TPC = remote.tpc( ...
            'imgVoltage',   ImgVoltage, ...
            'imgCurrent',   ImgCurrent, ...
            'pushVoltage',  PushVoltage, ...
            'pushCurrent',  PushCurrent, ...
            0);

        % add a bfoc acmo and built the sequence
        bfocId = 1;
        Sequence = usse.usse( 'acmo',Acmo, 'TPC',TPC, 'Popup',0 );
        Sequence = Sequence.setParam( 'Repeat',repetition, 'DropFrames',DropFrames );

        [Sequence, NbAcq] = Sequence.buildRemote;

        % Sequence.plot_diagram; xlim( [ 0 20e5 ] ); return

        %% build beamforming luts if needed
        BFStruct.Info.DebugMode = 1 ;
        disp('LUT computation');
        BFStruct = processing.lutBfoc( Sequence, bfocId, ImgInfo, NThreads );
%         if system.probe.Radius > 0
%             Data = processing.lutScanConvertImage( BFStruct, ImgInfo, 1024, 1024 ); % TODO: remove hardcoded size
%         end

        %% Setup RF buffer
        all_bufs = {};
        buf_size = Sequence.InfoStruct.mode(1).channelSize(1) * system.hardware.NbRxChan / 2 ;
        all_bufs{1} = int16( zeros( 1, buf_size, 'int16' ) );

        BFStruct.Info.DebugMode = 0 ;

        %   Initialize & load sequence remote
        Sequence = Sequence.initializeRemote( 'IPaddress',IPaddress, 'InitSequence',1 );
        % rmpath /home/labo/git_master/us/legHAL/src/../lib/libRemoteV9
        % addpath ../../rubi/remoteClients/matlab/remoteBuild/

        Sequence = Sequence.loadSequence();

        % Start sequence
        disp('start');
        Sequence = Sequence.startSequence();

        clear buf
        time_to_get_data = zeros( 1, NbAcq*repetition );
        time_to_loop = time_to_get_data;

        for k=1:NbAcq*repetition
            k
            tic_loop = tic;
            disp('Getting buffer');

            buffer_data = all_bufs{1};

            tic_get_data = tic;
            buf{k} = Sequence.getData('Realign',0, 'Timeout',120, 'Buffer',buffer_data, 'DataFormat',Sequence.getParam('DataFormat'), 'Debug',Debug );
            time_to_get_data(k) = toc(tic_get_data);

            % Stop sequence
            disp('reconstruction');

            tic_recon = tic;
            BFStruct = processing.reconBFoc( BFStruct, buf{k} );
            time_to_recon = toc(tic_recon);
    
            X = ImgInfo.XOrigin  + ( 0:(ImgInfo.nbTx * ImgInfo.RxLinesPerTx - 1) ) * ( system.probe.Pitch * ImgInfo.txPitch / ImgInfo.RxLinesPerTx ) ;
            Z = ImgInfo.Depth(1) + ( 0:(BFStruct.Lut(1).NPixLine - 1 ) ) * 1e-3 * common.constants.SoundSpeed / (RxFreq/4) * ImgInfo.PitchPix ;
            Z(1)
            Z(end)
            img =  abs( BFStruct.IQ ) ;

            time_to_loop(k) = toc(tic_loop);

            if Debug
                disp( [ 'Total time to loop : ' num2str(time_to_loop(k)) ] )
            end
        end
        Sequence = Sequence.stopSequence('Wait', 0);
        if DutyCycleLoop == 0
            img_noisefloor = img;
            noise_rf = buf;
        else
            img_brightness = img;
            brightness_rf = buf; 
        end
    end
 
    figure(11)
    imagesc(X,Z,20*log10(img_brightness(:,:,1)));
    title('aixplorer beamformed image')
    axis equal ;
    axis tight ;
    caxis ([72 135]);
    colormap gray ;
    
    size(img_brightness)

    % brightness
    for z = 1:size(img_brightness,1)
        brightness(z) = 20*log10(mean(img_brightness(z,:,1)));
    end

    % noise floor
    for z = 1:size(img_noisefloor,1)
        noiseFloor(z) = 20*log10(mean(img_noisefloor(z,:)));
        noiseStDev(z) = 20*log10(std(img_noisefloor(z,:)));
    end

end
