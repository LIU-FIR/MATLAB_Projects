%%
clear;clc;
fc = 1e6;
df = 1e5;
f0 = 0;
srate = 10e6;
T = 1e-3;
c = df/T;
t = 0:1/srate:10e-3;
% x_chirp = sin(2*pi*fc*t + 2*pi*c*t.^2/2);
% x_chirp = sin(2*pi*fc*t + 2*pi*(c*t.^2/2 + f0*t));
x_chirp = sin(2*pi*fc*t + 2*pi*(c*t.^2/2 + f0*t));
x_chirp = sin(2*pi*(c*t.^2/2 + f0*t));
%%
figure(1),clf 
% subplot(121)
plot(x_chirp,'k-','linew',2)
xlim([0 10000])
% subplot(122)
%%
figure(2), clf
pspectrum(x_chirp,srate,'spectrogram','FrequencyResolution',60e3)
%%

phi_t = 2*pi*(c*t.^2/2 + f0*t);
figure(3)
plot(t,phi_t,'b-','linew',1)
% xlim([0 .5e-3])