%%
clear;clc;

M = 512; % PFB Channels
R = 1; % Windows FIR Taps
N = M*R; % Overall N-FFT length

K = 16; % Zoom FFT-length
L = N*K; % length of the signal.
D = M*K; % 

n = 0:L-1; 
k = 0:M/2;

A = 2;
SNR_dB = -10;
SNR_LIN = db2mag(SNR_dB);
disp(['SNR in dB is: ',num2str(SNR_dB)])
disp(['SNR in linearity is: ',num2str(SNR_LIN)])
%

srate = 2000e6; % sampling rate
fv = 488.5e6; % signal frequency (Hz)
% fv = 455.8e6; % signal frequency (Hz)
% fv = 460e6; % signal frequency (Hz)
% fv = 488037109.375;


kvN = N*fv/srate; % frequency bin @ N-FFT
kvN_round = round(kvN);
disp(['N-pnt-FFT-frequency bin: ',num2str(kvN)])
disp(['Integral Npnt-FFT-frequency bin: ',num2str(kvN_round)])
% disp(['PFFT-frequency bin: ',num2str(kvm)])
% sig = A*exp(1i*2*pi*kv/N.*n);
sig = A*cos(2*pi*kvN/N*n); % carrier signal.
noise = (std(sig)/SNR_LIN)*randn(1,L); % zero-mean, SNR_LIN

x = sig + noise;

%%
% Compute M-FFT directly, with L/M frames.

frames_M_FFT = floor(L/M);

XM = fft(reshape(x,[M,frames_M_FFT]))/M;
XM_power = abs(XM).^2;

XM_power_avg = zeros(M,floor(frames_M_FFT/R));
for j = 1:floor(frames_M_FFT/R)
    XM_power_avg(:,j) = mean(XM_power(:,(j-1)*R+1:j*R),2);
end

idv_peak_XM = find(XM_power_avg(1:M/2+1)==max(XM_power_avg(1:M/2+1)));
kv_peak_XM = idv_peak_XM - 1;  % signal frequency-bin

XM_kvs = XM(idv_peak_XM,1:K); % select the signal channel.
Q = fft(XM_kvs)/K;
Q_power = abs(Q).^2;

idq_peak = find(Q_power==max(Q_power)); 
kq_peak = idq_peak - 1; 
if kq_peak > K/2
    kq_peak = kq_peak - K;
end

fr_M_FFT = (srate/M*(kv_peak_XM) + srate/M/K*kq_peak)/1e6;

%%
% Compute D=M*K FFT directly, with L/D frames;

frames_D_FFT = floor(L/D);
XD = fft(reshape(x,[D,frames_D_FFT]))/D;
XD_power = abs(XD).^2;
XD_power_avg = mean(XD_power,2);

idv_peak_XD = find(XD_power_avg(1:D/2+1)==max(XD_power_avg(1:D/2+1)));
kv_peak_XD = idv_peak_XD - 1;  % signal frequency-bin

fr_D_FFT = srate/D * kv_peak_XD/1e6;
disp(['The recovered signal-frequency via M-FFT: ',num2str(fr_M_FFT),' MHz'])
disp(['The recovered signal-frequency via D-FFT: ',num2str(fr_D_FFT),' MHz'])
%%
figure(1),clf
plot(0:D/2, pow2db(XD_power_avg(1:D/2+1)./max(XD_power_avg)),'k-','linew',2)
title([num2str(D),'-FFT power response'])
xlabel('FFT-bin')
ylabel('Normalized power: (dB)')
set(gca,'ylim',[-70 5])
grid on
%%
figure(2),clf
subplot(211)
plot(k,pow2db(XM_power_avg(1:M/2+1,1)./max(XM_power_avg(:,1))),'b-','linew',2)
title([num2str(M),'-Coarse FFT power response'])
xlabel('Coarse FFT-bin')
ylabel('Normalized power: (dB)')
grid on
set(gca,'ylim',[-20 5])
subplot(212)
plot(-K/2:K/2-1,fftshift(pow2db(Q_power)),'ko-','linew',2)
title([num2str(K),'-Zoom FFT power response'])
xlabel('Zoom FFT-bin: (zero-centered)')
ylabel('power: (dB)')
set(gca,'ylim',[-50 5])
grid on
% plot(k,X_power(1:M/2+1),'k-','linew',2)
%%
DPO = mod(kvN_round,R); % downsampled-phase-offset
os_ratio_str = '1';
os_ratio = str2num(os_ratio_str); % the oversampled ratio. os_ratio>1 means that outputs from PFFT's bins refresh faster.
slide_len = N/os_ratio; % when being loaded,the start position in orig.sequence increase by slide_len
frames_PFFT = floor((L-N)/slide_len)+1; % deducted from geometric relation

disp(['Downsampled-phase-offset(DPO): ',num2str(DPO)])
disp(['M-PFFT frames: ' num2str(frames_PFFT)])

rect_win = (rectwin(N))';
ham_win = (hamming(N))'; % hamming window function
sinc_win = fir1(N-1,1/M,rectwin(N)); % sinc window function
sinc_win_norm = sinc_win./max(sinc_win); % nomralized sinc window function
sys_win = sinc_win_norm.*ham_win; % composed system window function

Yp_frames = zeros(M,frames_PFFT);

for cnt = 1:frames_PFFT
    % for every sliding, compute the current head and tail of the
    % ready-to-pfft sequence.
    head_ptr = 1+(cnt-1)*slide_len;
    tail_ptr = 1+(cnt-1)*slide_len + N-1;
    
    x_reg = x(head_ptr:tail_ptr); % load the ready-to-pfft sequence in register.
%     Yp_frames(:,cnt) = pfft(x_reg,rect_win,M,R);
    Yp_frames(:,cnt) = pfft_n(x_reg,rect_win,M,R,DPO);
%     Yp_frames(:,cnt) = pfft_n(x_reg,sys_win,M,R,DPO);

end

%


Yp_power_avg = mean(abs(Yp_frames).^2,2);

% search for the signal frequency-bin M-index.
idv_peak = find(Yp_power_avg(1:M/2+1)==max(Yp_power_avg(1:M/2+1))); 
kv_peak = idv_peak - 1;  % signal frequency-bin



Yp_kvs = Yp_frames(idv_peak,1:K); % select the signal channel.

Z = fft(Yp_kvs)/K;
Z_power = abs(Z).^2;

% find 2nd-FFT's peak frequency.
idz_peak = find(Z_power==max(Z_power)); 
kz_peak = idz_peak - 1; 
if kz_peak > K/2
    kz_peak = kz_peak - K;
end

disp(['Peak PFFT-frequency bin: ',num2str(kv_peak)])
disp(['Peak 2nd-FFT frequency bin: ',num2str(kz_peak)])
fr = (srate/M*(kv_peak) + srate/N*DPO + srate/M/K*kz_peak)/1e6;
disp(['The original signal-frequency: ',num2str(fv/1e6),' MHz'])
disp(['The recovered signal-frequency via PFFT: ',num2str(fr),' MHz'])
disp(['The recovered signal-frequency via M-FFT: ',num2str(fr_M_FFT),' MHz'])
disp(['The recovered signal-frequency via D-FFT: ',num2str(fr_D_FFT),' MHz'])
%%
kz = -K/2:K/2-1;
figure(3),clf
subplot(211)
plot(k,pow2db(Yp_power_avg(1:M/2+1)./max(Yp_power_avg(1:M/2+1))),'b-','linew',2)
grid on
ylabel('Power (dB)')
set(gca,'ylim',[-20 0])
subplot(212)
plot(kz,fftshift(pow2db(Z_power./max(Z_power))),'ko-','linew',2)
grid on
ylabel('Power (dB)')

%%

