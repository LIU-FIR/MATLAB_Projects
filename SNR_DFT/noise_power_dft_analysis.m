%%
% 产生物理噪声信号，经过DFT后分析其显示功率谱密度和原功率谱密度的关系
% 1) 考虑矩形窗的(1/N)*DFT
clear;clc;
L = 1024;
Nfft = 1024;
M = 1000;
N = Nfft*M;
n= (0:N-1);

fs = 250e6;
fres = fs/Nfft;

% BW = fres;
% BW = enbw(hann(L))*fres;
% Pn_measured_dBm = 17; % 噪声功率测量值13 dBm

sigma = 0.05;% w/Hz

x = (fs)^0.5*sigma.*randn(1,N); % 在fs采样后的噪声序列，其功率为sigma^2*fs

x_arr = reshape(x,Nfft,M);
%%
% fcn_hanning = hann(L);
fcn_win = rectwin(L); % 可使用其他窗hann(L)
% fcn_hanning_arr = repmat(fcn_hanning,1,M);
fcn_win_arr = repmat(fcn_win,1,M);
x_win_arr = x_arr.*fcn_win_arr;

kfrex = 0:Nfft/2;
kfind = kfrex+1;
X_arr = fft(x_win_arr)/Nfft; % normalized fft
% PY_1 = abs(Y_arr(kfind,1)).^2;
% Y_mean = fft(x_arr)/Nfft; % normalized fft
Pxx_mean = mean(abs(X_arr).^2,2);

FNG = Nfft/2; FNG_dB = pow2db(FNG);
fres_dB = pow2db(fres); 

Pxx = Pxx_mean(kfind)*2;
Pxx_psd = Pxx/fres;
Pxx_dB = pow2db(Pxx);
Pxx_psd_dB = pow2db(Pxx_psd);
Pxx_psd_dB1 = Pxx_dB - fres_dB;
% noverlap = 0;
% [xpsd,f] = pwelch(x,L,noverlap,Nfft,fs);

%%
Pxx_esm = mean(Pxx(2:end-1));
sigma_esm = (Pxx_esm/2/fres)^0.5;
%%
figure(1)
subplot(2,1,1)
plot(kfrex*fres, Pxx_dB,'k-')
xlabel('Frequency (Hz)')
ylabel('PS (dB)')
subplot(2,1,2)
plot(kfrex*fres,Pxx_psd_dB,'b-')
hold on
plot(kfrex*fres,Pxx_psd_dB1,'r--')
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')

%%

% CPG_measure_db = pow2db(max(xpsd))-max(PY_PSD_dB);
% CPG_hanning = (sum(hann(L))/Nfft)^2;
% CPG_gain_dB_hanning = 10*log10(CPG_hanning);
% ICPG_hanning = sum(hann(L).^2)/Nfft;
% ICPG_gain_dB_hanning = 10*log10(ICPG_hanning);

% kc = find(xps==max(xps));
% xps_dbm = pow2db(xps/1e-3);
% xps_dbm(kc-1:kc+1) = pow2db(xps(kc-1:kc+1)/1e-3/ohm);
% 
% xps_mw = xps/1e-3;
% noise_xps_mean_dbm = pow2db(mean([xps_mw(2:kc-10);xps_mw(kc+10:end-1)]));
% sig_dbm = pow2db(xps(kc)/1e-3/ohm);
% figure(2)
% plot(f,xps_dbm,'b-')
% 
% disp([sig_dbm,noise_xps_mean_dbm])
