
% get data from old acmo
% Cano.RxDuration = A{1}.getParam('RxDuration');
% Cano.RxDelay = A{1}.getParam('RxDelay');
% Cano.RxFreq = A{1}.getParam('RxFreq');
% Cano.TwFreq = A{1}.getParam('TwFreq');

AixplorerIP    = 'local'; 

use_hadamard = 1;

Repeat_seq = 1;

Cano.RxDuration = 50;
Cano.RxDelay = 0;
Cano.TwFreq = 9;
Cano.RxFreq = 4*Cano.TwFreq;
Cano.NbHcycle = 2;

Cano.TGC = 400;

Cano.TxDelays = 0;

Cano.TxElemts = 1:system.probe.NbElemts;
Cano.RxElemts = 1:system.probe.NbElemts;

% Create canonical elusev
ELUSEV = elusev.hadamard( ...
... ELUSEV = elusev.canonical( ...
    'TrigIn', 0, ...
    'TrigOut', 0, ...
    'TrigAll', 0, ...
    'TxElemts', Cano.TxElemts, ...
    'TxDelays', Cano.TxDelays, ... % scalar: all, vector: size of TxElemts
    'TwFreq', Cano.TwFreq, ...
    'NbHcycle', Cano.NbHcycle, ...
    'DutyCycle', 1, ...
    'RxElemts', Cano.RxElemts, ... % -1: emitting element, 0: all, [...]>0 only given elements
    'RxFreq', Cano.RxFreq, ...
    'RxDelay', Cano.RxDelay, ...
    'RxDuration', Cano.RxDuration, ...
    'VGAInputFilter', 0, ...
    'RxBandwidth', 1, ...
    'FIRBandwidth', 100, ...
    'Pause', 5, ...
    'Repeat', 1, ...
    'PauseEnd', 8, ...
    0);

% ACMO for the sequence
ACMO = acmo.acmo( ...
    'elusev', ELUSEV, ...
    'ControlPts', Cano.TGC, ...
    0);

% USSE for the sequence
SEQ = usse.usse( ...
    'acmo', ACMO, ...
    ... 'Repeat', Repeat_seq, ...
    'TPC', remote.tpc('imgVoltage', 70, 'imgCurrent', 1), ...
    'Popup',0, ...
    0);

% Build remote structure
[Sequence, NbAcq] = SEQ.buildRemote();

% return

% Sequence execution

% Initialize & load sequence remote
Sequence = Sequence.initializeRemote( 'IPaddress', AixplorerIP, 'InitTGC',Cano.TGC, 'InitRxFreq',Cano.RxFreq );
Sequence = Sequence.loadSequence();

% disp('press a key !')
% pause

clear buf buf_128
for r = 1:Repeat_seq
    % Start sequence
    Sequence = Sequence.startSequence( 'Wait', 1 );

    % Retrieve data
    disp( 'getting data' )
    tic
    for k = 1 : NbAcq
        buf_128(r,k) = Sequence.getData( 'Realign', 2 ); %, 'RxChans', RxElemts );
    end
    toc

    % Stop sequence
    Sequence = Sequence.stopSequence( 'Wait', 1 );
    
    buf{r} = buf_128(r).RFdata(1:end-128,:,1:2:end) + buf_128(r).RFdata(1:end-128,:,2:2:end);
end


return

%% below is to check

%% hadamard
NbTx = length( TxElemts );
Hadamard_Matrix = hadamard( NbTx );
Hadamard_Matrix_inv = inv( Hadamard_Matrix );

data_hadamard = buf{1};
data_cano = zeros( size( data_hadamard ), 'single' );
for n = 1:NbTx
    n
    for n2 = 1:NbTx
        n2_coef = Hadamard_Matrix_inv( n2, n );
        data_cano(:,:,n) = data_cano(:,:,n) + single(data_hadamard(:,:,n2) ) * n2_coef; 
    end
end

%% plot acquisition
% data = data_cano ; fig_nb = 124124;
data = buf{1} ; fig_nb = 124126;
for n = 1:size(data,3)
    figure(fig_nb)
    imagesc( squeeze( data(:,:,n) ) )
    caxis( [-1 1] * 2e1 )
    
    pause(0.1)
end



%% below is to remove maybe
return


%%
r=1;
buf_tmp2 = buf_tmp(r).RFdata(1:end-128,:,1:2:end) + buf_tmp(r).RFdata(1:end-128,:,2:2:end);

buf = zeros(size(buf_tmp2));
for r = 1:Repeat_seq
    buf_tmp2 = buf_tmp(r).RFdata(1:end-128,:,1:2:end) + buf_tmp(r).RFdata(1:end-128,:,2:2:end);
    buf = buf + double(buf_tmp2);
end
buf = buf / Repeat_seq;

figure; imagesc( buf(:,:,8) )

return

%% imp
imp = zeros( [ size(buf.data,1)-128 system.probe.NbElemts ] );
for n = 1:length(TxElemts)
    imp(:,TxElemts(n)) = squeeze( sum( buf.data(1:end-128,TxElemts(n),[2*n-1 2*n] ), 3 ) );
end

figure; imagesc( imp ); caxis( [-1 1]*1e4 )
% figure;plot( imp(:,8) )

return

elec_HS_RX = [ 2 3 4 310  81 374   155 411  199 204 212  455 460 468 ] ;

%% BTE impedance


% using synthetic acquisition
NbRepeat = ELUSEV.getParam('Repeat');
NbTx = length( TxElemts );
RFsum = sum( buf.data(1:end-128, :, 1:NbRepeat*NbTx*2 ), 3 );

figure; imagesc( RFsum )
caxis([ -1 1 ] *1e4 )

return

%%
NbTx = length(TxElemts);
NbRx = length(RxElemts);
NbSmp = size( buf.data, 1 )-128;
NbRF = size( buf.data, 3 );
RFmean = zeros( [ NbSmp NbTx NbRx ], 'single' );

for n=1:NbTx
    n
    RFmean(:,:,n) = mean( buf.data(1:end-128,:,[(2*n-1):2*NbRx:NbRF (2*n):2*NbRx:NbRF]), 3 ) ;
end


figure; imagesc( sum( RFmean, 3 ) ); caxis( [-1 1]*1e3 )

return

%%
figure
for n=1:NbTx
    imagesc( RFmean(:,:,n) )
    caxis( [-1 1]*1e3 )
    title(num2str(n))
    pause %(0.01)
end

%%


[B,A] = butter(4,0.5/Sequence.InfoStruct.rx(1).Freq*2,'high');
RF_filt = filter(B,A,single(buf.data),[],3);

figure; 
imagesc(RF_filt(:,:,1))

% for n=1:128; imagesc( buf.data(:,:,n) ); caxis( [-1 1]*5e3); pause( 0.1 ); end


%%
load ~/buffer_hadamard.mat

TxElemts = 1:128;
NbTx = length( TxElemts );
Hadamard_Matrix = hadamard( NbTx );
Hadamard_Matrix_inv = inv( Hadamard_Matrix );

data2 = zeros( size( buffer.data ), 'single' );
for n = 1:NbTx
    n
    for n2 = 1:NbTx
        n2_coef = Hadamard_Matrix_inv( n2, n );
        data2(:,:,n) = data2(:,:,n) + single(buffer.data(:,:,n2) ) * n2_coef; 
    end
end

%%
load ~/buffer_canonic.mat
data1 = buffer.data;

figure( 762362 )
for n = 1:NbTx
    subplot( 2,1, 1)
    imagesc( data1(:,:,n) ); caxis( [-1 1]*1e3 *8);

    subplot( 2,1, 2)    
    imagesc( data2(:,:,n) );
    caxis( [-1 1]*3 * 128 * 8 )
    title( num2str(n) )
    pause( 0.1 )
end





