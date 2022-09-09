%%
% 产生物理信号，经过DFT后和Scoll给出的示例对应
clear;clc;
L = 1024;
Nfft = 1024;
M = 1000;
N = Nfft*M;
n= (0:N-1);

fs = 250e6;
fres = fs/Nfft;
k0 = 180;
w0 = 2*pi*k0/Nfft;


% BW = fres;
BW = enbw(hann(L))*fres;
ohm = 50; % 50 ohm
Psig_measured_dBm = 15; % 信号功率测量值11.9dBm
Pn_measured_dBm = 17; % 噪声功率测量值13 dBm

V0 = (10^(Psig_measured_dBm/10)*1e-3*ohm*2)^0.5;
sigma = (10^(Pn_measured_dBm/10)*1e-3/BW/2) ^0.5;


% SNR_BW_dB = 10*log10(V0^2/2/(sigma^2*BW));

s = V0*cos(w0*n); %产生周期信号, 幅度为电压值 (Volt)

wgnoise = (fs)^0.5*sigma.*randn(1,N); 
x = s + wgnoise;
x_arr = reshape(x,Nfft,M);
%%
fcn_hanning = hann(L);
fcn_boxcar = rectwin(L);
fcn_hanning_arr = repmat(fcn_hanning,1,M);
fcn_boxcar_arr = repmat(fcn_boxcar,1,M);
% x_win_arr = x_arr.*fcn_hanning_arr;
x_win_arr = x_arr.*fcn_hanning_arr;

kfrex = 0:Nfft/2;
kfind = kfrex+1;
Y_arr = fft(x_win_arr)/Nfft; % normalized fft
% PY_1 = abs(Y_arr(kfind,1)).^2;
% Y_mean = fft(x_arr)/Nfft; % normalized fft
PY_mean = mean(abs(Y_arr).^2,2);

PY_mean_SB = PY_mean*2; % 单边带single-band乘以2

fft_gain = 10*log10(Nfft/2);
enbw_hz = enbw(hann(L))*fres;

PY_PSD_hz = PY_mean_SB/enbw_hz;
PY_PSD_dB = pow2db(PY_PSD_hz); 
PY_PS_dB = pow2db(PY_PSD_hz*enbw_hz);

PY_PSD_hz_SD = PY_PSD_hz(kfind); % 仅保留单边
PY_PSD_dB_SD = PY_PSD_dB(kfind);
PY_PS_dB_SD = PY_PS_dB(kfind);
noverlap = 0;
[xpsd,f] = pwelch(x,hann(L),noverlap,Nfft,fs);
% [xpsd,f] = pwelch(x,L,noverlap,Nfft,fs);
%%
xps = enbw_hz.*xpsd;
CPG_measure_db = pow2db(max(xpsd))-max(PY_PSD_dB);

CPG_hanning = (sum(hann(L))/Nfft)^2;
CPG_gain_dB_hanning = 10*log10(CPG_hanning);
ICPG_hanning = sum(hann(L).^2)/Nfft;
ICPG_gain_dB_hanning = 10*log10(ICPG_hanning);

% pow2db(mean(xpsd(Nfft/4:end-1)))-pow2db(mean(PY_PSD_hz_SD(Nfft/4:end-1)))
% 
% pow2db(max(xps)/1e-3/ohm)
% pow2db(mean(xps(Nfft/4:end-1))/1e-3)


figure(1)
subplot(2,2,1)
plot(kfrex*fres, PY_PSD_dB_SD,'k-')
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
subplot(2,2,3)
plot(f,pow2db(xpsd),'b-')
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')

subplot(2,2,2)
plot(kfrex*fres, PY_PS_dB_SD,'k-')
xlabel('Frequency (Hz)')
ylabel('PS (dB)')
subplot(2,2,4)
plot(f,pow2db(xps),'b-')
xlabel('Frequency (Hz)')
ylabel('PS (dB)')

kc = find(xps==max(xps));
% xps_sig_dbm = pow2db(max(xps)/1e-3/ohm);
xps_dbm = pow2db(xps/1e-3);
xps_dbm(kc-1:kc+1) = pow2db(xps(kc-1:kc+1)/1e-3/ohm);

xps_mw = xps/1e-3;
noise_xps_mean_dbm = pow2db(mean([xps_mw(2:kc-10);xps_mw(kc+10:end-1)]));
sig_dbm = pow2db(xps(kc)/1e-3/ohm);
figure(2)
plot(f,xps_dbm,'b-')

disp([sig_dbm,noise_xps_mean_dbm])
% subplot(2,1,1)
% subplot(2,1,2)

% plot(kfrex, PY_PS_dB(kfind),'b-')

% noise_bw_gain_dB = 10*log10(fres*2);
% % PY_1_bw_dB = 10*log10(PY_1) + noise_bw_gain_dB;
% % PY_1_bw_dB(k0+1-1:k0+1+1) = 10*log10(PY_1(k0+1-1:k0+1+1));
% 
% 
% % PY_mean_bw_dB = 10*log10(PY_mean_SB)+ noise_bw_gain_dB;
% PY_mean_bw_dB = 10*log10(PY_mean_SB/1e-3)+noise_bw_gain_dB; % 标定为dBm
% % PY_mean_bw_dB(k0+1-1:k0+1+1) = 10*log10(PY_mean_SB(k0+1-1:k0+1+1));
% % PY_mean_bw_dB(k0+1-1:k0+1+1) = 10*log10(PY_mean_SB(k0+1-1:k0+1+1)/ohm/1e-3); % 标定为dBm
% 
% kc = k0+1; kspan = 1; klim = 20;
% noise_mean_power = mean([PY_mean_bw_dB(1:kc-klim);PY_mean_bw_dB(kc+klim:end)]);
% 
% ind1 = find(abs(PY_mean_bw_dB(1:end)-noise_mean_power)>2);
% PY_mean_bw_dB(ind1) = noise_mean_power;
% 
% PY_mean_bw_dB(kc-kspan:kc+kspan) = 10*log10(PY_mean_SB(kc-kspan:kc+kspan)/ohm/1e-3); % 中心3个点标定为dBm
% noise_bw_gain_dB1 = 10*log10(BW1);
% PY_1_bw_dB1 = 10*log10(PY_1) + noise_bw_gain_dB1;
% PY_1_bw_dB1(k0+1) = 10*log10(PY_1(k0+1));
%%



%%
% 用pwelch进行功率谱估计
clc;clear;
fs = 250e6;
% fc = 25e6;

M = 1000;
Nfft = 1024;
L = Nfft;
N = L*M;
fres = fs/Nfft;
enbw_hz = enbw(hanning(L))*fres;

BW = fres;


k0 = 180;
w0 = 2*pi*k0/Nfft;
n = 0:N-1;
ohm = 50; % 50 ohm
Psig_measured_dBm = 13; % 信号功率测量值11.9dBm
Pn_measured_dBm = 14; % 噪声功率测量值13 dBm

V0 = (10^(Psig_measured_dBm/10)*1e-3*ohm*2)^0.5;
sigma = (10^(Pn_measured_dBm/10)*1e-3/BW/2) ^0.5;

% s = V0*cos(w0*n); %产生周期信号, 幅度为电压值 (Volt)
x = V0*cos(w0*n) + sigma*(fs/2)^0.5*randn(size(n));

enbw_gain = 10*log10(enbw_hz);
bw_gain = 10*log10(fs/Nfft);
fft_gain = 10*log10(Nfft/2);
%%
noverlap = 0;
[xpsd,f] = pwelch(x,hanning(L),noverlap,Nfft,fs);
% [xps,f1] = pwelch(x,hanning(L),noverlap,Nfft,fs,'power');
xps = enbw_hz.*xpsd;

figure(1)
subplot(2,1,1)
plot(f,pow2db(xpsd))
xlabel('Frequency (Hz)')
ylabel('PSD (dB/Hz)')
subplot(2,1,2)
plot(f,pow2db(xps))
xlabel('Frequency (Hz)')
ylabel('PS (dB)')


%%
noise_sp = mean(xps(Nfft/2:end-1));
noise_dbm = pow2db(noise_sp/1e-3);
sig_dbm = pow2db(max(xps)/ohm/1e-3);
disp([sig_dbm,noise_dbm,Pn_measured_dBm-noise_dbm, enbw_gain-bw_gain]);
%%
sig_spd_dbm = pow2db(max(xpsd));
sig_sp_dbm = pow2db(max(xps));
disp(sig_sp_dbm - sig_spd_dbm)
%%

