% To test the symmetry of the pulse for THI on Aixplorer

clear all
close all
clc

addpath('/home/labo/git/us/oscilloscope/matlab')
addpath('/home/labo/git/us/oscilloscope/bin')
cd /home/labo/git/us/legHAL/
addPathLegHAL

vol = 65;
freq = 4.5;
duty_cycle = 1;
nb_half_cycle = 2;
save_dir = '/home/labo/git/us/legHAL/data';
file_name = Filenames('/home/labo/git/us/legHAL/data/');

for up = [false, true]

    % Initialize the oscilloscope
    trigg_level = 1;
    trigg_channel = 2;
    ip = '192.168.4.154';
    set_single_trigg_scope( ip, trigg_channel, trigg_level );
    % Test pulse polarity
    % up = false
    % Tramsit signal
    scripts.fire(vol, nb_half_cycle, freq, duty_cycle, up, 1, 256);

    % Read waveform from the oscilloscope
    read_channel = 1;
    TIME_OUT = 1000; %[ms]
    [ w, t ] = get_waveform( ip, read_channel, TIME_OUT );

    % Check data

    figure
    plot(t, w)
    save_name = file_name.generate('wfm_no_matching', {'V', 'pol'}, ...
                                [uint8(vol),  double(up)], 'full', true);

    save([save_name '.mat'], 't', 'w')
    pause(.4)

end

close all




