clear all ; clear classes ;
% cd /home/labo/git/us/legHAL
% addPathLegHAL
addpath( fullfile('..', 'tools', 'color_signal') )

% general parameters
modeRF = 0 ;
do_recon = 1;
NThreads = 4;
IPaddress = '192.168.3.125'; % IP address of Aixplorer

repetition = 1 ; % repeat sequence [1 inf]
% Voltage parameters
ImgVoltage  = 50;             % imaging voltage 0V]
ImgCurrent  = 0.1;            % security limit for imaging current [A]
PushVoltage = 10;             % imaging voltage [V]
PushCurrent = 5;              % security limit for imaging current [A]


% Image info
% fondamental
ImgInfo.TxPolarity     = [ 0 ];
ImgInfo.TxElemsPattern = [ 0 ];

ImgInfo.SteeringAngle  = 0; % steering angle [not implemented yet]
ImgInfo.XOrigin        = 0; % xmin [mm]
ImgInfo.FirstLinePos   = 86 % position of the first line in elt [>=1]
ImgInfo.txPitch        = 2; % txpitch [elements] (dist between tx in element pitch)
ImgInfo.nbTx           = 12 % txpitch [elements]
ImgInfo.RxLinesPerTx   = 4; % number of received line per tx (beamforming)
ImgInfo.PitchPix       = 0.5; % dz / lambda
ImgInfo.Depth(1)       = 40; % initial depth [mm]
ImgInfo.Depth(2)       = 30 + ImgInfo.Depth(1); % final depth [mm]
ImgInfo.RxApod(1)      = 0.4 % RxApodOne
ImgInfo.RxApod(2)      = 0.6; % RxApodZero
ImgInfo.EnsembleLength = 128;

ImgInfo.pulseInversion = 0;

% Emission parameters
TwFreq          = 3;      % transmit frequency [MHz]
NbHcycle        = 6;         % nb tx half cycle
ApodFct         = 'none';    % Apod function on the transmit aperture
TxFDRatio       = 3;         % F/D ratio (transmit)
DutyCycle       = 1;         % tx Duty cycle [0-1]
Focus           = 60 ;       % focal depth [mm]

% Color Parameters
Scale = 25; % cm/s
NumberSkippedFramesBeforeWallFilter = 1;
CutOffWallFilter = 0.03;
SmoothKernelSize = 0.4; %Size of the isotropic gaussian kernel to smooth R0 and R1 [mm]
CPIMinDisplay = 77; %dB
CPIMaxDisplay = 97; %dB

% ============================================================================ %
PRF    = 4*Scale*10/(common.constants.SoundSpeed/1000/TwFreq);         % Hz

% Receive parameters
RxBandwidth  = 2;          % 1 = 200%, 2 = 100%
RxFilterBandwidth = 50 ;
% TGC = [300 400 500 600 700 800];
TGC = [800 800 800 800 800 800];

%% converting ImgInfo parameters into acmo parameters
RxFreq       = 4 * TwFreq * (1+ImgInfo.pulseInversion); % MHz

NumSamples   = 128+ceil(1.4 * (2 * ImgInfo.Depth(2) - ImgInfo.Depth(1)) * 1e3 ...
    * RxFreq / (common.constants.SoundSpeed * 2^(RxBandwidth-1) * 128)) * 128;
SkipSamples  = floor(0. * ImgInfo.Depth(1) * RxFreq * 1e3 ...
    / common.constants.SoundSpeed);

RxDuration = NumSamples / RxFreq * 2^(RxBandwidth-1);
RxDelay    = SkipSamples / RxFreq;

%% built sequence

Acmo = acmo.color(...
    'TwFreq',       TwFreq, ...
    'NbHcycle',     NbHcycle, ...
    'SteerAngle',   ImgInfo.SteeringAngle, ...
    'Focus',        Focus, ...
    'PRF',          PRF, ...
    'DutyCycle',    1, ...
    'XOrigin',      ImgInfo.XOrigin, ...
    'FirstLinePos', ImgInfo.FirstLinePos,...
    'NbFocLines',   ImgInfo.nbTx,...
    'TxLinePitch',  ImgInfo.txPitch,...
    'FDRatio',      TxFDRatio, ...
    'ApodFct',      'none', ...
    'RxFreq',       RxFreq, ...
    'RxDuration',   RxDuration, ...
    'RxDelay',      RxDelay, ...
    'RxBandwidth',  RxBandwidth, ...
    'FIRBandwidth', RxFilterBandwidth, ...
    'PulseInv',     ImgInfo.pulseInversion, ...
    'TxPolarity',   ImgInfo.TxPolarity, ...
    'TxElemsPattern', ImgInfo.TxElemsPattern, ...
    'ControlPts',   TGC, ...
    'EnsembleLength', ImgInfo.EnsembleLength, ...
    'NbLocalBuffer',2, ...
    'NbHostBuffer', 2, ...
    0);

TPC = remote.tpc( ...
    'imgVoltage',   ImgVoltage, ...
    'imgCurrent',   ImgCurrent, ...
    0);

% adds a bfoc acmo and builds the sequence
colorId = 1;
Sequence = usse.usse( 'acmo',Acmo, 'TPC',TPC, 'Popup', 0 );
if ~strcmp(system.probe.Name, 'SC6-1')
    warning('Please select the SC6-1 probe, then run "clear classes" and run the script again');
    Sequence = Sequence.selectProbe();
end
Sequence = Sequence.setParam( 'Repeat',repetition, 'DropFrames', 0 );

[Sequence NbAcq Acmo] = Sequence.buildRemote;

ImgInfo.NumberLinesPerBox = Acmo{1}.getParam('NumberLinesPerBox');
ImgInfo.nbTx =  Acmo{1}.getParam('NbFocLines');

% built beamforming luts if needed
if ~modeRF
    disp('LUT computation');
    BFStruct = processing.lutColor( Sequence, colorId, ImgInfo, NThreads );

%     if system.probe.Radius > 0
%         Data = processing.lutScanConvertImage( BFStruct, ImgInfo, 1024, 1024 ); % TODO: remove hardcoded size
%     end
end

% clear buffer if necessary
clear buffer;ImgInfo

BFStruct.Info.DebugMode = 0 ;

% return

%   Initialize & load sequence remote
Sequence = Sequence.initializeRemote( 'IPaddress',IPaddress, 'InitSequence',0 );
Sequence = Sequence.loadSequence();

%%
% Start sequence
disp('start');
Sequence = Sequence.startSequence();

for k=1:repetition
    disp('Getting buffer');

    if modeRF
        buffer = Sequence.getData('Realign', 2);
    else
        buffer = Sequence.getData('Realign', 0);
        
        if do_recon
            Nz = size(BFStruct.IQ,1);
            IQMatrix = zeros(Nz,ImgInfo.nbTx*ImgInfo.RxLinesPerTx,ImgInfo.EnsembleLength,repetition);

            disp('reconstruction');

            BFStruct = processing.reconColor( BFStruct, buffer );
            NumberBoxes = round(ImgInfo.nbTx/ImgInfo.NumberLinesPerBox);
            JIndexBoxZero = repmat(0:(ImgInfo.NumberLinesPerBox-1),[1 ImgInfo.EnsembleLength]);
            JIndex = [];
            for idxBox=0:(NumberBoxes-1)
                BoxOffsetLineNumber = idxBox * ImgInfo.NumberLinesPerBox;
                JIndex = [JIndex (JIndexBoxZero + BoxOffsetLineNumber)];
            end;%idxBox

            Nz = size(BFStruct.IQ,1);
            for idx3 = 0:(ImgInfo.nbTx*ImgInfo.EnsembleLength-1)
                j = JIndex(idx3+1);
                idxBox = floor(double(idx3)/double(ImgInfo.EnsembleLength*ImgInfo.NumberLinesPerBox));
                firingnumbercurrentbox = (idx3 - idxBox * double(ImgInfo.EnsembleLength) * double(ImgInfo.NumberLinesPerBox));
                PRIidx = floor(((    double(idx3) - double(idxBox) * double(ImgInfo.EnsembleLength) * double(ImgInfo.NumberLinesPerBox))/double(ImgInfo.NumberLinesPerBox)));
                IQMatrix(:,(1+j*ImgInfo.RxLinesPerTx):(ImgInfo.RxLinesPerTx+j*ImgInfo.RxLinesPerTx),PRIidx+1,k) = BFStruct.IQ(:,(idx3*ImgInfo.RxLinesPerTx+1):(idx3*ImgInfo.RxLinesPerTx+ImgInfo.RxLinesPerTx));
            end;


%             for idx=1:ImgInfo.EnsembleLength
% 
%                 if system.probe.Radius > 0
%                     Data2 = processing.ScanConvertImage( Data, BFStruct.IQ );
%                     img = Data2.img ;
%                     X = Data2.XAxis;
%                     Z = Data2.ZAxis;
%                 else
%                     X = ImgInfo.XOrigin  + ( 0:(ImgInfo.nbTx * ImgInfo.RxLinesPerTx - 1) ) * ( system.probe.Pitch * ImgInfo.txPitch / ImgInfo.RxLinesPerTx ) ;
%                     Z = ImgInfo.Depth(1) + ( 0:(BFStruct.Lut(1).NPixLine - 1 ) ) * 1e-3 * common.constants.SoundSpeed / (RxFreq/4) * ImgInfo.PitchPix ;
%                     img = 20*log10( abs( IQMatrix(:,:,idx,k) ) );
%                 end

%                 figure(10);
%                 imagesc(X,Z,img);
%                 axis equal ;
%                 axis tight ;
% %                 caxis( [0 40]+110 )
%                 colormap gray ;
%                 title(['idx = ',num2str(idx),'; k = ',num2str(k)])
% 
%                 drawnow ;
%                 pause(.1);
%             end
        end
    end;
    % clear buffer;
end
Sequence = Sequence.stopSequence('Wait', 1);

%%
figure(52352)
imagesc( 20*log10( abs( squeeze(IQMatrix(:,:,1)) ) ) )
colormap gray
colorbar
caxis( [ 100 160 ] )

% return

%%
if modeRF
    nb_box = ImgInfo.nbTx / ImgInfo.NumberLinesPerBox;
    for n_box = 1:nb_box
        for n_el = 1:ImgInfo.EnsembleLength
            for n_line = 1:ImgInfo.NumberLinesPerBox
                n = (n_box-1)  * ImgInfo.NumberLinesPerBox * ImgInfo.EnsembleLength ...
                    + (n_el-1) * ImgInfo.NumberLinesPerBox ...
                    + n_line ;
                figure(124141); imagesc( buffer.data(:,:,n) );
                title( [ 'box ' num2str(n_box) ', n el ' num2str(n_el) ' , line ' num2str(n_line) ] )
                pause
            end
        end
    end
    return
end

%Wall Filter computation
% projection initialization method

[B,A] = butter(3, CutOffWallFilter*2, 'high');
order=3;
el=ImgInfo.EnsembleLength-NumberSkippedFramesBeforeWallFilter;
g=B((order+1):-1:2)-B(1)*A((order+1):-1:2);
g = g';
F=eye(order-1);
F=horzcat(zeros(order-1,1),F);
F=vertcat(F,-A(order+1:-1:2));
%
B_init=zeros(el,order);
for n=0:el-1
    B_init(n+1,:)=g'*(F^n);
end
d=B(1);
q=zeros(order-1,1);
q=vertcat(q,[1]);
c=eye(el)*d;
cnt=0;
%
for j=1:length(c)
    for i=1:length(c)
        if (i>j)
            c(i,j)=g'*((F^cnt)*q);
            cnt=cnt+1;
        end
    end
    cnt=0;
end
H = (eye(el)-(B_init*inv((B_init'*B_init))*B_init'))*c;

%% spectrum
addpath /home/labo/git/us/tools/color_signal/
Fe = single(Acmo{1}.getParam('PRF')); ... 
    .../ single(Acmo{1}.getParam('NumberLinesPerBox')) ;

IQSpectra_normalized(IQMatrix, Fe, 1, 100);

return

%%
for k=1:repetition
    IQ=IQMatrix(:,:,(1+NumberSkippedFramesBeforeWallFilter):ImgInfo.EnsembleLength,k);
    [d1,d2,d3] = size(IQ);
    Y = zeros(d1,d2,d3);
    for i1 = 1:d1
        for i2 = 1:d2
            x=  IQ(i1,i2,:);
            x = reshape(x,[el,1,1]);
            y = H*x;
            Y(i1,i2,:) = reshape(y,[1,1,el]);
        end;
    end;
    t0=1;
    t1=el;
    
    R1 = mean(conj(Y(:,:,t0:(t1-1))).*Y(:,:,(t0+1):t1),3);
    
    
    dx = system.probe.Pitch * ImgInfo.txPitch / ImgInfo.RxLinesPerTx;
    dz = 1e-3 * common.constants.SoundSpeed / (RxFreq/4) * ImgInfo.PitchPix;
    Nj = ceil(3*SmoothKernelSize/dx);
    Ni = ceil(3*SmoothKernelSize/dz);
    [I,J]= ndgrid(-Ni:Ni,-Nj:Nj);
    ZKernel = dz * I;
    XKernel = dx * J;
    KernelMatrix = exp((-XKernel.^2-ZKernel.^2)/2/SmoothKernelSize/SmoothKernelSize);
    KernelMatrix = KernelMatrix / sum(KernelMatrix(:));
    
    
    clear CFI BMode CPI
    
    FilteredR1 = filter2(KernelMatrix,R1);
    CFI(:,:,k) = angle(FilteredR1)/2/pi;
    BMode(:,:,k) = mean(abs(IQMatrix(:,:,:,k)),3);
    R0 = 10*log10(mean(abs(Y(:,:,:)).^2,3));
    CPI(:,:,k) = filter2(KernelMatrix,R0);
    
    II = (0:31)'/31;O = zeros(32,1);
    VelocityColormap = [[flipud(II),O,O];[O,O,II]];
    figure(3);
    imagesc(X,Z,CFI(:,:,k),[-0.5 0.5]);
    colormap(VelocityColormap);
    title(['CFI; FrameIdx = ',num2str(k)]);
    axis equal;axis tight;
    
    bmin=77;
    bmax=97;
    figure(2);
    colormap(hot);
    imagesc(X,Z,CPI(:,:,k),[CPIMinDisplay CPIMaxDisplay]);
    title(['CPI; FrameIdx = ',num2str(k)]);
    
    axis equal;axis tight;
end;

AverageCFI = mean(CFI,3);
AverageCPI = mean(CPI,3);
AverageBMode = 20*log10(mean(BMode,3));

figure(11);
colormap(gray);
imagesc(X,Z,AverageBMode);
title('AverageBMode');


figure(13);
imagesc(X,Z,AverageCFI,[-0.5 0.5]);
colormap(VelocityColormap);
title('Average velocity');
axis equal;axis tight;


figure(12);
colormap(hot);
imagesc(X,Z,AverageCPI,[CPIMinDisplay CPIMaxDisplay]);
title('Average CPI');


return


