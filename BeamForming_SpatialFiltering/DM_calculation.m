%%
% compute pulse broadening width

clear; clc;

f_H = 346e6;
f_L = 314e6;
fc = (f_H + f_L)/2/1e6; % MHz
Nch = 1024;
B = (f_H - f_L)/1e6; % MHz
df = B/Nch;
DM = 100; % dispersion measurement

D = 4.15e6; % MHz^2 PC^(-1)cm^3 ms
%%
dW = 2*D*df/fc^3 * DM;

disp([' width of the pulse with the intrinsic width W:', num2str(dW)])

%%
clear;clc;
%%
srate=1000; f=[0 250];
t=0:1/srate:2; n=length(t);
tmp = linspace(f(1),f(2),n);
% chirpTS = sin(2*pi.*linspace(f(1),f(2)*mean(f)/f(2),n).*t);
chirpTS = sin(2*pi.*tmp(end:-1:(n-1)/2+1).*t(1:(n+1)/2));

% t = 0:1/1e3:1;
% y = chirp(t,0,1,250);
%%
pspectrum(chirpTS,srate,'spectrogram',...
    'FrequencyResolution',25,'OverlapPercent',99,'Leakage',0.85)

