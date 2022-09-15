%%
% Compute S21 response 
%

clear;clc;

Qi = 348000;
Qc = 43000;
df = -5e3; % Hz
% df = 0;
Qr = 1/(1/Qi+1/Qc);
% Qr = 17000;
% Qc = 20000;
fr = 5e9; % resonance frequency
fc = 1.3e9; % carrier frequency
fs = 4e9;
tn = 0:1/fs:10e-6;
xc = exp(1i*2*pi*fc*tn);

BW = 1e6;
fres = 1e3;
f = fr-BW/2:fres:fr+BW/2;
fk = f-fr;
%%
S_21 = 1-Qr/Qc*(1-1i*2*Qc*df/fr)./(1+1i*2*Qr*(f-fr)./fr);
%
dly = 0.75e-2;
theta = dly/3e8*fr*2*pi;
loss = .5^0.5;

S_21_m = exp(1i*theta).*S_21*loss;
S_21_m_dB = mag2db(abs(S_21_m));
S_21_m_pha = angle(S_21_m);

S_21_dB = mag2db(abs(S_21));
S_21_pha = angle(S_21);
%%
figure(1),clf
subplot(211),hold on
plot(fk,S_21_dB,'b-','linew',2)
plot(fk,S_21_m_dB,'m-','linew',2)
subplot(212),hold on
plot(fk,S_21_pha,'b-','linew',2)
plot(fk,S_21_m_pha,'m-','linew',2)

%%
% y = S_21(101).*xc;


ctr_idx = floor(length(f)/2)+1;
figure(2),clf
hold on
plot(real(S_21),imag(S_21),'b.')
plot(real(S_21_m),imag(S_21_m),'m.')
set(gca,'xlim',[-1.2 1.2],'ylim',[-1.2 1.2])
grid on, hold on, axis square
plot(get(gca,'xlim'),[0 0],'k','linew',1)
plot([0 0],get(gca,'ylim'),'k','linew',1)
plot(real(S_21(ctr_idx)),imag(S_21(ctr_idx)),'rx','markersize',8,'markerfacecolor','r','linew',2)
plot(real(S_21_m(ctr_idx)),imag(S_21_m(ctr_idx)),'kx','markersize',8,'markerfacecolor','r','linew',2)
%%

