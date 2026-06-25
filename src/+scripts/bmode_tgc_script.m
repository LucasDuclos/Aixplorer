 

% paths of created files
folder = '/home/labo/git/us/legHAL/data';
if exist(folder,'dir')==0
    mkdir(folder)
end

ip = '192.168.4.0';

voltage = 50;
duty_cycle = 0.66; % [0 1]

demodFreq = 7.5;           % MHz
filterBw = 80;         % % 
speedOfSound  = 1540;     % m/s
lambda = speedOfSound / demodFreq / 1000;    % mm
depth = 250*lambda;              % mm
%vec_TGC = [100 240 380 520 660 800 920];

p.txFreq = 7.5;
p.nb_halfCycles = 2;

p.txLinePitch = 2;
p.nbTXChannels = 128; 
p.zOrigin = 1.15;
p.focusMm = 2*depth/3;            % mm
p.steeringAngle = 0;       % °
bfg.nbLines = 4;
bfg.pixelPitch = 0.5;

step = 100;
vectgc = 0:step:900;
tgc_cpt = 1;
for idxtgc = vectgc 
    tgc_val = idxtgc;
    vec_TGC = ones(1,7)*tgc_val;
    [Z, noiseFloor, noiseStDev, brightness, mean_roi, mean_roi2 ] = scripts.bmode_tgc(   ip,...
                                                                    voltage,...
                                                                    duty_cycle,...
                                                                    demodFreq,...
                                                                    filterBw,...
                                                                    depth,...
                                                                    vec_TGC,...
                                                                    p,...
                                                                    bfg    );

   
   mean_brightness(tgc_cpt) = mean_roi;
   mean_noisefloor(tgc_cpt) = mean_roi2;
   tgc_cpt = tgc_cpt + 1;
end
      
figure(18)
plot(vectgc, mean_noisefloor)
title('noisefloor vs tgc, 23 mm <= z <= 28 mm')
xlabel('Tgc')
ylabel('Mean Intensity (dB)')
figure
plot(vectgc, mean_brightness)
title('brightness vs tgc, 23 mm <= z <= 28 mm')
xlabel('Tgc')
ylabel('Mean Intensity (dB)')

% figure
% plot(brightness,Z)
% hold on
% plot(noiseFloor, Z, 'r')
% set(gca,'YDir','reverse')
% title('Aixplorer bmode brightness and noise floor');
% legend('brightness', 'noise floor', 'Location', 'southeast')
% grid on
% 
% figure
% plot(brightness-noiseFloor,Z)
% hold on
% set(gca,'YDir','reverse')
% title('SNR versus depth');
% grid on
% 
% fileName = [folder '/dc' num2str(duty_cycle) '_dmod' num2str(demodFreq) '_fBw' num2str(filterBw) '_tgc' num2str(length(vec_TGC)) '.mat'];
% save(fileName,'noiseFloor', 'noiseStDev', 'brightness');