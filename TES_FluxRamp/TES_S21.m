%%
% Compute S21 response 
%

clear;clc;

Qi = 348000;
Qc = 43000;
% df = -5e3; % Hz
df = 0;
Qr = 1/(1/Qi+1/Qc);
fr = 5e9; % resonance frequency
fc = 1.3e9; % carrier frequency
fs = 4e9;
tn = 0:1/fs:10e-6;
xc = exp(1i*2*pi*fc*tn);

BW = 2e6;
fres = 10e3;
f = fr-BW/2:fres:fr+BW/2;
fk = f-fr;
%%
S_21 = 1-Qr/Qc*(1-1i*2*Qc*df/fr)./(1+1i*2*Qr*(f-fr)./fr);
%
S_21_dB = mag2db(abs(S_21));
S_21_pha = angle(S_21);
%%
figure(1),clf
subplot(211)
plot(fk,S_21_dB,'k-','linew',1)
subplot(212)
plot(fk,S_21_pha,'b-','linew',1)

%%
y = S_21(101).*xc;
figure(2),clf
plot(real(y),imag(y),'b.')