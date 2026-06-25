% Probe Charac main script
clear all; close all; clc;

Charac.AixplorerIP = '192.168.3.107';
% Select probe
probe = 'SL15_4';
% Output file name;

matching_value = inputdlg('matching value (uH)');
matching_value = cell2mat(matching_value);
VoltageVect = [10 20 30];
TGCVect = [600 600 600];
file_prefix = 'Test';
SaveData_Path ='/home/labo/TX_RX_Characterization/KoloTest1'; 

common.constants.SoundSpeed = 1540;
if (strcmp(probe, 'SL15_4') == 1)
    Charac.NPiezos        = 256;
    Charac.FrequencyRange = [4 15];
    Charac.TGC            = 600;
    Charac.RxDelay = 5;
    Charac.RxDuration = 40;
elseif (strcmp(probe, 'NV1D_256') == 1)
    Charac.NPiezos        = 256;
    Charac.FrequencyRange = [4 15];
    Charac.TGC            = 600;
    Charac.RxDelay = 1;
    Charac.RxDuration =50;
elseif (strcmp(probe, 'L10_2') == 1)
    Charac.NPiezos        = 192;
    Charac.FrequencyRange = [2 10];
    Charac.TGC            = 600;
    Charac.RxDelay = 5;
    Charac.RxDuration =40;
elseif (strcmp(probe, 'LV16_5') == 1)
    Charac.NPiezos        = 192;
    Charac.FrequencyRange = [5 16];
    Charac.TGC            = 600;
    Charac.RxDelay = 5;
    Charac.RxDuration =50;
elseif  (strcmp(probe, 'P5_1X') == 1)
    Charac.NPiezos        = 80;
    Charac.TGC            = 500;
    Charac.FrequencyRange = [1 6];
    Charac.RxDelay = 5;
    Charac.RxDuration =10;
elseif  (strcmp(probe, 'LH20_6') == 1)
    Charac.NPiezos        = 192;
    Charac.TGC            = 600;
    Charac.FrequencyRange = [4 20];
    Charac.RxDelay = 0;
    Charac.RxDuration =40;
elseif  (strcmp(probe, 'SC6_1') == 1)
    Charac.NPiezos        = 192;
    Charac.TGC            = 600;
    Charac.FrequencyRange = [1 6];
    Charac.RxDelay = 60;
    Charac.RxDuration =50;
elseif  (strcmp(probe, 'E12_3') == 1)
    Charac.NPiezos        = 192;
    Charac.TGC            = 600;
    Charac.FrequencyRange = [3 12];
    Charac.RxDelay = 10;
    Charac.RxDuration =50;
elseif (strcmp(probe, 'VL20_6') == 1)
    Charac.NPiezos        = 256;
    Charac.FrequencyRange = [4 15];
    Charac.TGC            = 600;
    Charac.RxDelay = 0;
    Charac.RxDuration =40;
elseif  (strcmp(probe, 'Oldelft_TEE_64') == 1)
    Charac.NPiezos        = 129;
    Charac.TGC            = 600;
    Charac.FrequencyRange = [3 12];
    Charac.RxDelay = 0;
    Charac.RxDuration =50;
else
    disp([probe ' probe not supported yet']);
end

Charac.NbAcq  = 1;
Charac.Voltage = 10;
Charac.RxFreq = 30;   % MHz unit
% Charac.RxDelay = 60;
% Charac.RxDuration =50;
Charac.Pause = 200;
Charac.PauseEnd = 2000;
Charac.Nfft = 4096;
Charac.AllChannReceive = 0;

% % %%%%% Define waveform %%%%%/home/labo/TX_RX_Characterization/VL20_6
waveform = [1 1 1 1 1 0 0 0];
% waveform = [0 0 0 0 0 0 0 0];
for nk = 1:Charac.NPiezos*2;
Charac.CharactElements(nk) = ceil(nk/2);
Charac.Waveform(nk, 1:length(waveform)) = (-1)^nk*waveform(:)';
end
[perChannelData infoStruct] = scripts.AcqProbeCharac(Charac);
%%% refine time window


SkipSamples = double(infoStruct(1).event(1).skipSamples);
nbSamplesDepth = double(infoStruct(1).event(1).numSamples);
SequenceRxFreq = double(infoStruct(1).rx(1).Freq);

DepthAxis = ((([1:nbSamplesDepth] +SkipSamples)./(SequenceRxFreq*1e6))*1480/2)*1e3;
figure(1);
% imagesc(1:Charac.NPiezos*2, DepthAxis, perChannelData);
% ylabel('mm')
% xlabel('element')
imagesc(perChannelData);
caxis( [-1 1]*1e2 );

%%
[X,Y] = ginput(2);
% 
% X = [187.6; 188.5];
% Y = [70; 807];

% update RxDelay and RxDuration
Charac.RxDelay = (infoStruct.event(1).skipSamples+min(Y))/Charac.RxFreq;
Charac.RxDuration = abs(Y(2)-Y(1))/Charac.RxFreq;
Charac.RxFreq = 60;
Charac.NbAcq = 10;

for measIndex = 1:length(VoltageVect)
    Charac.Voltage = VoltageVect(measIndex);
    Charac.TGC = TGCVect(measIndex);
    [perChannelData infoStruct] = scripts.AcqProbeCharac(Charac);
    
    perChannelDataAvg = mean(double(perChannelData),3);
    perChannelDataStd = std(double(perChannelData),1,3);
    
    %%%% Plot pulse averaged over all acq of the central element
    figure;
    plot(0:(size(perChannelData,1)-1),perChannelDataAvg(:,end/2));
    title(['Average Pulse middle of the probe / Voltage: ', num2str(Charac.Voltage) ,' V'])
    
    %%%%%  Imagesc Pulses averaged over all acq and standart deviation aver
    %%%%%  all acquisition
    figure;
    subplot(2,1,1);
    imagesc(perChannelDataAvg); title(['Average Pulse / Voltage: ', num2str(Charac.Voltage) ,' V']);
    colorbar
    subplot(2,1,2);
    imagesc(perChannelDataStd); title(['Standard deviation over NbAcq / Voltage: ', num2str(Charac.Voltage) ,' V']);
    colorbar
    
    %%%%% Compute Spactrum 
    spectrums = fft(double(perChannelData),Charac.Nfft,1)/Charac.Nfft;
    spectrumsAvg = abs(mean(spectrums,3));
    spectrumsStd = std(spectrums,1,3);
    
    spectrumsAvgAcrossArray = mean(spectrumsAvg,2);
    spectrumsStdAcrossArray = std(spectrumsAvg,1,2);
    
    % per channel bandwidth, fc, cuttoff
    
    idxfftMin = round(Charac.FrequencyRange(1)*Charac.Nfft/Charac.RxFreq);
    idxfftMax = round(Charac.FrequencyRange(2)*Charac.Nfft/Charac.RxFreq);
    
    %%%%% Compute Bandwidth info (BW,central freq, High and Low cutoff) of each element (averaged over all Acquisition) 
    for k=1:length(Charac.CharactElements)
        [Maxfft,idxMax] = max(spectrumsAvg(idxfftMin:idxfftMax,k));
        fMax(k) = (idxMax+idxfftMin-1)/Charac.Nfft*Charac.RxFreq;
        aboveM6dB = find(spectrumsAvg(idxfftMin:idxfftMax,k)>=0.5*Maxfft);
        lowCutOff6dB(k) = (min(aboveM6dB)+idxfftMin-1)*Charac.RxFreq/Charac.Nfft;
        highCutOff6dB(k) = (max(aboveM6dB)+idxfftMin-1)*Charac.RxFreq/Charac.Nfft;
    end
    %%%%% Statistics of BW 
    % Mean and std of central freq
    meanFc = mean(lowCutOff6dB+highCutOff6dB)/2;
    stdFc = std((lowCutOff6dB+highCutOff6dB)/2);
    
    % Mean and std of low cutoff
    meanLowCutOff6dB = mean(lowCutOff6dB);
    stdLowCutOff6dB = std(lowCutOff6dB);
    
    % Mean and std of high cutoff
    meanHighCutOff6dB = mean(highCutOff6dB);
    stdHighCutOff6dB = std(highCutOff6dB);
    
    % Mean and std of Bandwidth
    meanBW = mean(2*(highCutOff6dB-lowCutOff6dB)./(highCutOff6dB+lowCutOff6dB));
    stdBW = std(2*(highCutOff6dB-lowCutOff6dB)./(highCutOff6dB+lowCutOff6dB));
    
    %Peak Peak Voltage center of probe
    CentralPiezo = round(Charac.NPiezos / 2);
    CentralRFAvg = perChannelDataAvg(:,CentralPiezo);
    PeakPeakRcv = max(CentralRFAvg) - min(CentralRFAvg);
    
    %-20 dB pulse duration
    CentralAvgEnveloppe = abs(hilbert(CentralRFAvg));
    CentralAvgEnveloppedB = 20*log10(CentralAvgEnveloppe);
    MaxEnvdB = max(CentralAvgEnveloppedB);
    idxMinus20dB = processing.FindZero(CentralAvgEnveloppedB - MaxEnvdB + 20);
    if (length(idxMinus20dB)>1)
        PulseDuration20dB_us = (max(idxMinus20dB) - min(idxMinus20dB)) / Charac.RxFreq;
    else
        PulseDuration20dB_us = 0.0;
    end;
    
    
    %%%%% Save statistics 
    stat.meanFc = meanFc;
    stat.stdFc = stdFc;
    stat.meanLowCutOff6dB = meanLowCutOff6dB;
    stat.stdLowCutOff6dB = stdLowCutOff6dB;
    stat.meanHighCutOff6dB = meanHighCutOff6dB;
    stat.stdLowCutOff6dB = stdLowCutOff6dB;
    stat.meanBW = meanBW;
    stat.stdBW = stdBW;
    stat.PeakPeakRcv = PeakPeakRcv;
    stat.PulseDuration20dB_us = PulseDuration20dB_us;
    
    % tof over the array
    channelRef = round(size(perChannelData,2)/2);
    for acq=1:size(perChannelData,3)
        for chann=1:size(perChannelData,2)
            Rxy = xcorr(perChannelData(:,chann,acq),perChannelData(:,channelRef,acq));
            [M,I] = max(abs(Rxy));
            tof(chann) = I-0.5*(Rxy(I+1)-Rxy(I-1))/(Rxy(I+1)+Rxy(I-1)-2*Rxy(I));
        end
        p = polyfit(Charac.CharactElements,tof,1);
        p2(acq,:) = p;
        residualTof(:,acq) = tof- polyval(p,Charac.CharactElements);
    end
    %
    misAligmentSlopeAvg = mean(1e3*p2(:,1)/Charac.RxFreq);
    misAligmentSlopeStd = std(1e3*p2(:,1)/Charac.RxFreq);
    %
    % probeMotionStd =  std((p2(2:end,2)-p2(1,2))/Charac.RxFreq);
    
    %%% Plot Spectra %%%%%
    figure;
    plot((0:(Charac.Nfft-1))/Charac.Nfft*Charac.RxFreq,20*log10(abs(spectrumsAvgAcrossArray)),'Linewidth',2);
    hold on;
    plot((0:(Charac.Nfft-1))/Charac.Nfft*Charac.RxFreq,20*log10(abs(spectrumsAvgAcrossArray)-3*spectrumsStdAcrossArray),'k--');
    plot((0:(Charac.Nfft-1))/Charac.Nfft*Charac.RxFreq,20*log10(abs(spectrumsAvgAcrossArray)+3*spectrumsStdAcrossArray),'k--');
    hold off;
    xlim(Charac.FrequencyRange);
    title(['Average spectrum over all acquisition and all elements / Voltage: ', num2str(Charac.Voltage) ...
        ,' V     fc: ', num2str(meanFc), ' MHz / LC: ', num2str(meanLowCutOff6dB), 'MHz / HC: ', ...
        num2str(meanHighCutOff6dB), ' MHz / BW: ', num2str(meanBW*100), ' %']);
    xlabel('Frequency MHz');
    
    
    %%%%% Plot Time of flight %%%%%
    figure;
    plot(Charac.CharactElements,mean(1e3*residualTof/Charac.RxFreq,2), '*');
    title(['time of flight / Voltage: ', num2str(Charac.Voltage) ,' V']);
    xlabel('channel');
    ylabel('time of flight (ns)');
    
    voltage = num2str(Charac.Voltage);
    %%%%% Display Bandwidth statistics info
    disp(['Voltage (V): ', num2str(Charac.Voltage)]);
    disp(['TGC: ', num2str(Charac.TGC)]);
    disp(['fc6dB (MHz): ' num2str(meanFc)]);
    disp(['BW6dB (BW/fc): ' num2str(meanBW)]);
    disp(['LCO6dB (MHz): ' num2str(meanLowCutOff6dB)]);
    disp(['HCO6dB (MHz): ' num2str(meanHighCutOff6dB)]);
    disp(['PeakPeakRcv: ' num2str(stat.PeakPeakRcv)]);
    disp(['PulseDuration20dB_us: ' num2str(stat.PulseDuration20dB_us)]);
    disp(['PeakPeakRcv / Voltage = ' num2str(stat.PeakPeakRcv / Charac.Voltage)]);
    
    
    cc = clock;
    data_time = [ num2str(cc(1)) '_' num2str(cc(2)) '_' num2str(cc(3)) '_' num2str(cc(4)) '_' num2str(cc(5)) ] ;
    
    if matching_value == 0
        inductUnit = '';
        matching_value = '';
    else
        inductUnit = 'uH_';
    end
    %%%%% Save Data
    save( [ SaveData_Path,'/', file_prefix '_' voltage 'V_ShortPulse_' num2str(matching_value) inductUnit,data_time '.mat' ], 'perChannelData','infoStruct' , 'Charac', 'stat');
end
