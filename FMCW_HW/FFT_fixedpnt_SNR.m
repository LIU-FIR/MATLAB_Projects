%%
clear;clc;

N = 8192; % FFT-length

D = 1; % FFT frames
L = N*D; % length of the signal.

n = (0:L-1)'; % for dsp.FFT

OutputWidth = 24;
OutputFractionWidth = 22;
ProductDataType = 'Custom';
AccumulatorDataType = 'Custom';
% ProductWidth = 18;
% ProductFractionWidth = 16;
% AccumulatorWidth = 18;
% AccumulatorFractionWidth = 16;

ProductWidth = 32;
ProductFractionWidth = 30;
AccumulatorWidth = 32;
AccumulatorFractionWidth = 30;


F = fimath('OverflowAction', 'Saturate', 'RoundingMethod',...
    'Floor', 'ProductMode','KeepMSB');
ft = dsp.FFT('FFTLengthSource','Property','FFTLength',N,'Normalize',true,...
    'OutputDataType','Custom','CustomOutputDataType',numerictype([],OutputWidth,OutputFractionWidth),...
    'ProductDataType','Custom','CustomProductDataType',numerictype([],18,16),...
    'AccumulatorDataType','Custom','CustomProductDataType',numerictype([],18,16),...
    'OverflowAction','Saturate');

FSV = 0.25; % Full-Scale-Volt 250mV, Vpp=500mV
A = 0.070; % signal's amplitude
SNR_dB = 1;
SNR_LIN = db2mag(SNR_dB);
% disp(['SNR in dB is: ',num2str(SNR_dB)])
% disp(['SNR in linearity is: ',num2str(SNR_LIN)])

srate = 1000e6; % sampling rate
fv = 250e6; % signal frequency (Hz)

kvN = N*fv/srate; % frequency bin @ N-FFT
idvN = kvN+1;
% disp(['N-pnt-FFT-frequency bin: ',num2str(kvN)])

sig = A*cos(2*pi*fv/srate*n); % carrier signal.

ntimes = 100; % 100 loops for simulation
SNR_dbl = zeros(ntimes,1);
SNR_fp = zeros(ntimes,1);
SNR_loss = zeros(ntimes,1);
%%
disp('Loop begins now')
for nt = 1:ntimes
    noise = (std(sig)/SNR_LIN)*randn(L,1); % zero-mean, std = std(sig)/SNR_LIN
    x = sig + noise;
    xs = x./FSV;
    %
    % figure(1),clf
    % subplot(211)
    % plot(x,'b-','linew',1)
    % subplot(212)
    % plot(xs,'k-','linew',1)  % Scale

    x_fp = fi(xs,1,12,11,F);
    % x_fp_dbl = double(x_fp);

    y_fp = ft(x_fp); % fixed-point-precision FFT
    % Y = abs(y).^2;
    y = fft(double(x_fp))/N; % double-precision FFT

    %
    Y_fp = abs(y_fp);
    Y = abs(y);

    %
    % figure(2),clf
    % subplot(211)
    % plot(1:L/2,pow2db(double(Y_fp(2:L/2+1))),'k-','linew',1)
    % subplot(212)
    % plot(1:L/2,pow2db(Y(2:L/2+1)),'b-','linew',1)
    %

    noise_fp_mean = mean([Y_fp(idvN-2:idvN-1).^2;Y_fp(idvN+1:idvN+2).^2]);
    noise_mean = mean([Y(idvN-2:idvN-1).^2;Y(idvN+1:idvN+2).^2]);
    SNR_fp(nt) = pow2db(double(Y_fp(idvN).^2/noise_fp_mean));
    SNR_dbl(nt) = pow2db(Y(idvN).^2/noise_mean);

    SNR_loss(nt) = SNR_dbl(nt) - SNR_fp(nt);

    % disp(['SNR is: ',num2str(SNR)])
    % disp(['SNR of fixed-point-math is: ',num2str(SNR_fp)])
    % disp(['SNR loss is: ',num2str(SNR_loss)])    
end
disp('end here')
%%
% release(ft);
figure(3),clf
plot(SNR_loss,'s','linew',1,'markerfacecol','r')
xlim([0 101])
% ylim([-0.1 0.1])
grid on
xlabel('Test times')
ylabel('dB')
title(['SNR loss: Product(',[num2str(ProductWidth),',',...
    num2str(ProductFractionWidth)],')',' Accumulator(',[num2str(AccumulatorWidth),',',...
    num2str(AccumulatorFractionWidth)],')',' Output(',[num2str(OutputWidth),',',...
    num2str(OutputFractionWidth)],')'])
%%
figure(1),clf
subplot(211)
plot(x,'b-','linew',1)
grid on
title('Data samples')
xlabel('Samples')
ylabel('Volt')
ylim([-FSV FSV])
subplot(212)
plot(xs,'m-','linew',1)  % Scale
grid on
title('Scaled data samples')
xlabel('Samples')
ylabel('Scaled data')
%

figure(2),clf
subplot(211)
plot(1:L/2,pow2db(Y(2:L/2+1).^2),'k-','linew',1)
title('Double-precision FFT power spectrum')
xlabel('FFT bins')
ylabel('dB')
grid on
subplot(212)
plot(1:L/2,pow2db(double(Y_fp(2:L/2+1)).^2),'r-','linew',1)
title('Fixed-point-precision FFT power spectrum')
xlabel('FFT bins')
ylabel('dB')
grid on
%%
