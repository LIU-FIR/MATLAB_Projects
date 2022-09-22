%%
% Compute S21 response 
%

clear;clc;
%%
Qi = 348000;
Qc = 43000;
% Qi = 41000;
% Qc = 129000;
df1 = 10e3; % Hz % asymmetry, resonance frequency deviation.
df = 0;
Qr = 1/(1/Qi+1/Qc);
% Qr = 17000;
% Qc = 20000;
fr = 5e9; % original resonance frequency
% fc = 1.3e9; % carrier frequency
% eta_offset = -0.05 + 0.1i;

% fs = 4e9;
% tn = 0:1/fs:10e-6;
% xc = exp(1i*2*pi*fc*tn);

BW = 300e3;
fstep = 5e2;
f = fr-BW:fstep:fr+BW;
fk = f-fr;
%
S =1-Qr/Qc*1./(1+1i*2*Qr*(f-fr)./fr); % Ideal S21
S_ms = 1-Qr/Qc*(1-1i*2*Qc*df1/fr)./(1+1i*2*Qr*(f-fr)./fr); % S21 measurement.
% S_21 = 1-Qr/Qc*(1-1i*2*Qc*df1/fr)./(1+1i*2*Qr*(f-fr)./fr);
%
dly = 1.2e-2;
theta = -dly/3e8*fr*2*pi;
loss = .75^0.5;

% S_21_m = (1-Qr/Qc*(1-1i*2*Qc*df1/fr)./(1+1i*2*Qr*(f-fr)./fr)).*loss.*exp(1i*theta);
% S_21_m_cal = (1-Qr/Qc*(1-1i*2*Qc*df1/fr)./(1+1i*2*Qr*(f-fr)./fr));
% S_21_m = exp(1i*theta).*S_21*loss;
S_ms_dB = mag2db(abs(S_ms));
S_ms_pha = angle(S_ms);
idx_ftone = find(abs(S)==min(abs(S))); % ftone is set at the original resonance frequency.
idx_fres = find(abs(S_ms)==min(abs(S_ms))); % shifted resonance frequency.

phi_cal = -angle(S_ms(idx_fres));
F_cal = exp(1i*phi_cal);

S_ms = F_cal*S_ms; % calibration to quadrature axis.
S_dB = mag2db(abs(S));
S_pha = angle(S);

%
Eta = (fk(idx_fres+1)-fk(idx_fres-1))./(S_ms(idx_fres+1)-S_ms(idx_fres-1));
disp(['Eta = ',num2str(Eta)]);
Eta1 = (fk(idx_ftone+1)-fk(idx_ftone-1))./(S_ms(idx_ftone+1)-S_ms(idx_ftone-1));
Eta2 = (fk(idx_fres+1)-fk(idx_fres-1))./imag(S_ms(idx_fres+1)-S_ms(idx_fres-1));

S_E_cal = Eta*S_ms;%/max(abs(Eta));
S_E_cal1 = Eta1*S_ms;%/max(abs(Eta));
S_E_cal2 = Eta2*S_ms;%/max(abs(Eta));

S_E_cal_dB = mag2db(abs(S_E_cal));
S_E_cal_pha = angle(S_E_cal);


Df_est = real(Eta*S_ms(idx_ftone));
Df_est1= real(Eta1*S_ms(idx_ftone));
Df_est2 = imag(Eta2*S_ms(idx_ftone));
disp([Df_est df1 Df_est1 Df_est2])

%%
figure(1),clf
% subplot(211),hold on
hold on
% plot(fk,S_dB,'k-','linew',2)
yyaxis left
plot(fk,abs(S_ms),'b-','linew',2)
% plot(fk,real(S_E_cal),'r-','linew',2)
% plot(fk,imag(S_E),'r-','linew',2)
% plot(fk,real(S_E_cal),'ms','linew',2)
plot([fk(idx_fres) fk(idx_fres)],[abs(S_ms(idx_fres)) abs(S_ms(idx_fres))],'ko','markersize',6,'markerfacecolor','k','linew',2)
xlabel('Hz')
ylabel('Magnitude: (a.u.)')
yyaxis right
plot(fk,S_ms_pha,'r-','linew',2)
ylabel('Phase: (rad)')
yrange = get(gca,'ylim');
xrange = get(gca,'xlim');
plot([fk(idx_fres) fk(idx_fres)],[yrange(1) yrange(2)],'k--','linew',2)
plot([fk(idx_fres) fk(idx_fres)],[S_ms_pha(idx_fres) S_ms_pha(idx_fres)],'ko','markersize',6,'markerfacecolor','k','linew',2)

% plot([xrange(1) xrange(2)],[0 0],'k-','linew',2)
grid on
% subplot(212),hold on

% plot(fk,S_E_cal_pha,'m-','linew',2)
% plot(fk,S_E_cal_pha-S_E_pha,'b-','linew',2)

%%
% y = S_21(101).*xc;



% ctr_idx = floor(length(f)/2)+1;
figure(2),clf
hold on
plot(real(S_ms),imag(S_ms),'b.')
grid on, hold on, axis square
set(gca,'xlim',[0 1.2],'ylim',[-.6 .6])
plot(get(gca,'xlim'),[0 0],'k','linew',1)
plot([0 0],get(gca,'ylim'),'k','linew',1)
plot(real(S_ms(idx_ftone)),imag(S_ms(idx_ftone)),'kx','markersize',8,'markerfacecolor','r','linew',2)
plot(real(S_ms(idx_fres)),imag(S_ms(idx_fres)),'rx','markersize',10,'markerfacecolor','r','linew',2)

plot([real(S_ms(idx_fres+10)) real(S_ms(idx_fres+10))],[imag(S_ms(idx_fres+10)) imag(S_ms(idx_fres+10))],'ms','markersize',8,'linew',2)
plot([real(S_ms(idx_fres-10)) real(S_ms(idx_fres-10))],[imag(S_ms(idx_fres-10)) imag(S_ms(idx_fres-10))],'ms','markersize',8,'linew',2)
plot([real(S_ms(idx_fres-10)) real(S_ms(idx_fres+10))],[imag(S_ms(idx_fres-10)) imag(S_ms(idx_fres+10))],'m-','linew',2)

figure(3),clf
hold on
% xlim([-4000 4000])
% ylim([-4000 4000])
% plot(real(S),imag(S),'k.')
% plot(real(S_E),imag(S_E),'b.')
plot(real(S_E_cal),imag(S_E_cal),'m.')
% plot(real(S_E_cal1),imag(S_E_cal1),'b.')
plot(real(S_E_cal2),imag(S_E_cal2),'b.')
% set(gca,'xlim',[-1.2 1.2],'ylim',[-1.2 1.2])
grid on, hold on, axis square
plot(get(gca,'xlim'),[0 0],'k','linew',1)
plot([0 0],get(gca,'ylim'),'k','linew',1)
% plot(real(S(ctr_idx)),imag(S(ctr_idx)),'rx','markersize',8,'markerfacecolor','r','linew',2)
% plot(real(S_E(ctr_idx)),imag(S_E(ctr_idx)),'kx','markersize',8,'markerfacecolor','r','linew',2)
plot(real(S_E_cal(idx_ftone)),imag(S_E_cal(idx_ftone)),'kx','markersize',8,'markerfacecolor','r','linew',2)
plot(real(S_E_cal(idx_fres)),imag(S_E_cal(idx_fres)),'rx','markersize',8,'markerfacecolor','r','linew',2)

% plot(real(S_E_cal1(idx_fr)),imag(S_E_cal1(idx_fr)),'kx','markersize',8,'markerfacecolor','r','linew',2)
% plot(real(S_E_cal1(idx_ftone)),imag(S_E_cal1(idx_ftone)),'rx','markersize',8,'markerfacecolor','r','linew',2)

plot(real(S_E_cal2(idx_ftone)),imag(S_E_cal2(idx_ftone)),'kx','markersize',8,'markerfacecolor','r','linew',2)
plot(real(S_E_cal2(idx_fres)),imag(S_E_cal2(idx_fres)),'rx','markersize',8,'markerfacecolor','r','linew',2)


%%

