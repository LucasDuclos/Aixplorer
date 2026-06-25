
% paths of created files
folder = '/home/labo/git/us/legHAL/data';
if exist(folder,'dir')==0
    mkdir(folder)
end

ip = '192.168.4.0';

voltage = 60;
duty_cycle = 0.66; % [0 1]

demodFreq = 7.5;           % MHz
filterBw = 80;         % %
speedOfSound  = 1540;     % m/s
lambda = speedOfSound / demodFreq / 1000;    % mm
depth = 250*lambda;              % mm

%vec_TGC = [405 473 525 562 641 734 814 930 930 930 930 930 930 930 930 930];
vec_TGC = [0 930];
p.txFreq = 7.5;
p.nb_halfCycles = 2;

p.txLinePitch = 2;
p.nbTXChannels = 128; 
p.zOrigin = 1.15;
p.focusMm = 2*depth/3;            % mm
p.steeringAngle = 0;       % °
bfg.nbLines = 4;
bfg.pixelPitch = 0.5;

[Z, noiseFloor, noiseStDev, brightness, noise_rf, brightness_rf  ] = scripts.bmode_snr(   ip,...
                                                                voltage,...
                                                                duty_cycle,...
                                                                demodFreq,...
                                                                filterBw,...
                                                                depth,...
                                                                vec_TGC,...
                                                                p,...
                                                                bfg    );
                                    
figure
plot(brightness,Z)
hold on
plot(noiseFloor, Z, 'r')
set(gca,'YDir','reverse')
%title({'Aixplorer bmode brightness and noise floor for' 'TGC = [405 473 525 562 641 734 814 930 930 930 930 930 930 930 930 930]'}, 'FontSize', 9);
title('Aixplorer bmode brightness and noise floor for a TGC ramp');
legend('brightness', 'noise floor', 'Location', 'southeast')
xlabel('Intensity (dB)')
ylabel('depth (mm)')
grid on

figure
plot(brightness-noiseFloor,Z)
hold on
set(gca,'YDir','reverse')
%title({'Aixplorer bmode SNR for'  'TGC = [405 473 525 562 641 734 814 930 930 930 930 930 930 930 930 930]'}, 'FontSize', 9);
title('Aixplorer bmode SNR for a TGC ramp');
xlabel('Intensity (dB)')
ylabel('depth (mm)')
grid on


fileName = [folder '/dc' num2str(duty_cycle) '_dmod' num2str(demodFreq) '_fBw' num2str(filterBw) '_tgc' num2str(length(vec_TGC)) '.mat'];
save(fileName,'noiseFloor', 'noiseStDev', 'brightness','noise_rf','brightness_rf');

%% deltaT avec burst
% [valmax, Imax] = max(brightness);
% [valmin, Imin] = min(brightness);
% 
% deltaz =  depth(Imax)-depth(Imin);
% deltaT = deltaz*eventDuration/depth;
