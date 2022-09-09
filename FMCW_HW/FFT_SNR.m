%%
clear;clc;

N = 8192; % FFT-length

D = 16; % FFT frames
L = N*D; % length of the signal.

n = 0:L-1; 
k = 0:N/2;

A = 2;
SNR_dB = -10;
SNR_LIN = db2mag(SNR_dB);
disp(['SNR in dB is: ',num2str(SNR_dB)])
disp(['SNR in linearity is: ',num2str(SNR_LIN)])
%

srate = 1000e6; % sampling rate
fv = 250e6; % signal frequency (Hz)
% fv = 455.8e6; % signal frequency (Hz)
% fv = 460e6; % signal frequency (Hz)
% fv = 488037109.375;


kvN = N*fv/srate; % frequency bin @ N-FFT
kvN_round = round(kvN);
disp(['N-pnt-FFT-frequency bin: ',num2str(kvN)])
disp(['Integral Npnt-FFT-frequency bin: ',num2str(kvN_round)])
% disp(['PFFT-frequency bin: ',num2str(kvm)])
% sig = A*exp(1i*2*pi*kv/N.*n);
sig = A*cos(2*pi*fv/srate*n); % carrier signal.
noise = (std(sig)/SNR_LIN)*randn(1,L); % zero-mean, std = std(sig)/SNR_LIN

x = sig + noise;

%%
% Compute N-FFT directly, with L/N frames;

frames_N_FFT = floor(L/N);
XN = fft(reshape(x,[N,frames_N_FFT]))/N;

XN_power = abs(XN).^2;
XN_power_avg = mean(XN_power,2);

idv_peak_XN = find(XN_power_avg(1:N/2+1)==max(XN_power_avg(1:N/2+1)));
kv_peak_XN = idv_peak_XN - 1;  % signal frequency-bin

fr_N_FFT = srate/N * kv_peak_XN/1e6;
disp(['The recovered signal-frequency via N-FFT: ',num2str(fr_N_FFT),' MHz'])
%%
figure(1),clf
subplot(211)
plot(0:N/2, pow2db(XN_power_avg(1:N/2+1)./max(XN_power_avg)),'k-','linew',2)
title([num2str(N),'-FFT average power response'])
xlabel('FFT-bin')
ylabel('Normalized power: (dB)')
set(gca,'ylim',[-70 5])
grid on
subplot(212)
plot(0:N/2, pow2db(XN_power(1:N/2+1,1)./max(XN_power(1:N/2+1,1))),'b-','linew',2)
title([num2str(N),'-FFT power response'])
xlabel('FFT-bin')
ylabel('Normalized power: (dB)')
set(gca,'ylim',[-70 5])
grid on
%% 
% Different fft-length 
clear;clc;
L = 16 * 8192;
srate = 1000e6; % sampling rate
fv = 250e6; % signal frequency (Hz)
A = 2;
SNR_dB = -10;
SNR_LIN = db2mag(SNR_dB);
n = 0:L-1; 
sig = A*cos(2*pi*fv/srate*n); % carrier signal.
noise = (std(sig)/SNR_LIN)*randn(1,L); % zero-mean, std = std(sig)/SNR_LIN

x = sig + noise;
% Nps = [512 1024 2048 4096 8192];
Nps = [512 8192];
%%
for k = 1:2
    N = Nps(k);
    frames_N_FFT = floor(L/N);
    kvN = N*fv/srate; % frequency bin @ N-FFT
%     kvN_round = round(kvN);
    disp(['N-pnt-FFT-frequency bin: ',num2str(kvN)])    
    XN = fft(reshape(x,[N,frames_N_FFT]))/N;

    XN_power = abs(XN).^2;
    XN_power_avg = mean(XN_power,2);

    idv_peak_XN = find(XN_power_avg(1:N/2+1)==max(XN_power_avg(1:N/2+1)));
    kv_peak_XN = idv_peak_XN - 1;  % signal frequency-bin

    fr_N_FFT = srate/N * kv_peak_XN/1e6;
    disp(['The recovered signal-frequency via N-FFT: ',num2str(fr_N_FFT),' MHz'])
    
    figure(k),clf
    subplot(211)
    plot(0:N/2, pow2db(XN_power_avg(1:N/2+1)./max(XN_power_avg)),'k-','linew',1)
    title([num2str(N),'-FFT average power response'])
    xlabel('FFT-bin')
    ylabel('Normalized power: (dB)')
    set(gca,'ylim',[-70 5])
    grid on
    subplot(212)
    plot(0:N/2, pow2db(XN_power(1:N/2+1,1)./max(XN_power(1:N/2+1,1))),'b-','linew',1)
    title([num2str(N),'-FFT power response'])
    xlabel('FFT-bin')
    ylabel('Normalized power: (dB)')
    set(gca,'ylim',[-70 5])
    grid on
end
