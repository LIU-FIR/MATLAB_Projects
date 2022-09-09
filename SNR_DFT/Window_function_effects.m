%%
% 计算coherent power gain 和 incoherent power gain
clear;clc;
L = 801;
Nfft = 2048;
fs = 10e3;
fres = fs/Nfft;

fcn_hanning = hann(L);
noise = 5*randn(1,L);

CG = sum(fcn_hanning)/L;
ENBW = L*sum(fcn_hanning.^2)/sum(fcn_hanning).^2;

CPG = 10*log10(CG^2);

disp([CG ENBW CPG])
%%
% 计算Fs=10 kHz,Nfft=2048, 501点长度的Hanning window的ENBW
clear;clc;

Nfft = 2048;
% L = 501;
L = Nfft;
fs = 10e3;
fres = fs/Nfft;

fcn_hanning = hann(L);
fcn_rect = rectwin(L);

enbw_gain_dB = 10*log10(enbw(hann(L)));

CPG_hanning = (sum(hann(L))/Nfft)^2;
ICPG_hanning = sum(hann(L).^2)/Nfft;

CPG_gain_dB_hanning = 10*log10(CPG_hanning);
ICPG_gain_dB_hanning = 10*log10(ICPG_hanning);
PG_gain_dB = CPG_gain_dB_hanning - ICPG_gain_dB_hanning;

disp([CPG_hanning,ICPG_hanning])
