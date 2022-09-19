%%
% Compute S21 response 
%

clear;clc;

Qi = 348000;
Qc = 43000;
% Qi = 41000;
% Qc = 129000;
df1 = 12e3; % Hz % asymmetry
df = 0;
Qr = 1/(1/Qi+1/Qc);
% Qr = 17000;
% Qc = 20000;
fr = 5e9; % resonance frequency
fc = 1.3e9; % carrier frequency
% eta_offset = -0.05 + 0.1i;

% fs = 4e9;
% tn = 0:1/fs:10e-6;
% xc = exp(1i*2*pi*fc*tn);

BW = 300e3;
fres = 2e3;
f = fr-BW:fres:fr+BW;
fk = f-fr;
%%
S =1-Qr/Qc*1./(1+1i*2*Qr*(f-fr)./fr); 
S_E = 1-Qr/Qc*(1-1i*2*Qc*df1/fr)./(1+1i*2*Qr*(f-fr)./fr);
% S_21 = 1-Qr/Qc*(1-1i*2*Qc*df1/fr)./(1+1i*2*Qr*(f-fr)./fr);
%
dly = 1.2e-2;
theta = -dly/3e8*fr*2*pi;
loss = .75^0.5;

% S_21_m = (1-Qr/Qc*(1-1i*2*Qc*df1/fr)./(1+1i*2*Qr*(f-fr)./fr)).*loss.*exp(1i*theta);
% S_21_m_cal = (1-Qr/Qc*(1-1i*2*Qc*df1/fr)./(1+1i*2*Qr*(f-fr)./fr));
% S_21_m = exp(1i*theta).*S_21*loss;
Eta = (fk(ind1+1)-fk(ind1-1))./(S_E(ind1+1)-S_E(ind1-1));
S_E_cal = Eta*S_E/max(abs(Eta));

S_E_dB = mag2db(abs(S_E));
S_E_pha = angle(S_E);

S_E_cal_dB = mag2db(abs(S_E_cal));
S_E_cal_pha = angle(S_E_cal);

S_dB = mag2db(abs(S));
S_pha = angle(S);

ind0 = find(abs(S)==min(abs(S)));
ind1 = find(S_E_dB==min(S_E_dB));
%%
figure(1),clf
subplot(211),hold on
% plot(fk,S_dB,'k-','linew',2)
plot(fk,S_E_dB,'b-','linew',2)
plot(fk,S_E_cal_dB,'ms','linew',2)
yrange = get(gca,'ylim');
plot([fk(ind0) fk(ind0)],[yrange(1) yrange(2)],'r-','linew',2)
subplot(212),hold on
plot(fk,S_E_pha,'k-','linew',2)
plot(fk,S_E_cal_pha,'m-','linew',2)
% plot(fk,S_E_cal_pha-S_E_pha,'b-','linew',2)

%%
% y = S_21(101).*xc;



ctr_idx = floor(length(f)/2)+1;
figure(2),clf
hold on
% plot(real(S),imag(S),'k.')
plot(real(S_E),imag(S_E),'b.')
plot(real(S_E_cal),imag(S_E_cal),'m.')
set(gca,'xlim',[-1.2 1.2],'ylim',[-1.2 1.2])
grid on, hold on, axis square
plot(get(gca,'xlim'),[0 0],'k','linew',1)
plot([0 0],get(gca,'ylim'),'k','linew',1)
plot(real(S(ctr_idx)),imag(S(ctr_idx)),'rx','markersize',8,'markerfacecolor','r','linew',2)
plot(real(S_E(ctr_idx)),imag(S_E(ctr_idx)),'kx','markersize',8,'markerfacecolor','r','linew',2)
plot(real(S_E_cal(ctr_idx)),imag(S_E_cal(ctr_idx)),'kx','markersize',8,'markerfacecolor','r','linew',2)

%
Df_est = real(Eta*S_E(ind0));
disp([Df_est df1])
