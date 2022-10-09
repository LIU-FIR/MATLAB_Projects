%%
% Compute S21 response 
%

clear;clc;
%%
Qi = 348000;
df1 = 10e3; % Hz % asymmetry, resonance frequency deviation.
df = 0;
% Qr = 1/(1/Qi+1/Qc); % relationship between Qr and Qi,Qc.
Qr = 17000;
Qc = 20000;
fres = 5e9; % original resonance frequency


BW = 100e3; % Yu SLAC Microresonator RF (SMuRF) Electronics, P8, CMB umux resonators bandwidth
fstep = 5e2; % computational frequency interval/step.
f = fres-BW/2*10:fstep:fres+BW/2*10; % computational frequency range, fres is central frequency. 
fzc = f-fres; % move fres to zero frequency,
%
S =1-Qr/Qc*1./(1+1i*2*Qr*(f-fres)./fres); % Ideal S21, fres at real axis.
% S_ms = 1-Qr/Qc*(1-1i*2*Qc*df1/fr)./(1+1i*2*Qr*(f-fr)./fr); % S21 measurement,asymmetry.
S_ms = S.*exp(1i*pi/8); % uncalibrated S21 measurement，fres is not located at real axis.


S_dB = mag2db(abs(S));
S_pha = angle(S);
S_ms_dB = mag2db(abs(S_ms));
S_ms_pha = angle(S_ms);

idx_offset = 10;

idx_fres = find(abs(S_ms)==min(abs(S_ms))); % "shifted" resonance frequency.
idx_ftone = idx_fres + idx_offset; %% assume ftone stays unchanged.

%
Eta = (fzc(idx_fres+1)-fzc(idx_fres-1))./(S_ms(idx_fres+1)-S_ms(idx_fres-1));
% Eta1 = (fk(idx_fres+1)-fk(idx_fres-1))./(S_ms(idx_fres-1)-S_ms(idx_fres+1));
disp(['Eta = ',num2str(Eta)]);
% disp(['Eta1 = ',num2str(Eta1)]);

S_ms_cal = Eta.*S_ms*1i;% calibrated S21 to vertical axis, Eta is a rotated vector.
S_ms_cal_v = Eta.*S_ms; % calibrated S21 to horizontal axis, Eta is a rotated vector.

Df_est = imag(S_ms_cal(idx_ftone));
disp([Df_est fzc(idx_ftone)])
%
%%
figure(1),clf
% subplot(211),hold on
hold on
yyaxis left
% plot(fk,abs(S),'k-','linew',2)
plot(fzc,abs(S_ms),'b-','linew',2)

plot([fzc(idx_fres) fzc(idx_fres)],[abs(S_ms(idx_fres)) abs(S_ms(idx_fres))],'ko','markersize',6,'markerfacecolor','k','linew',2)
xlabel('Hz')
ylabel('Magnitude: (a.u.)')
yyaxis right
plot(fzc,S_ms_pha,'r-','linew',2)
ylabel('Phase: (rad)')
yrange = get(gca,'ylim');
xrange = get(gca,'xlim');
plot([fzc(idx_fres) fzc(idx_fres)],[yrange(1) yrange(2)],'k--','linew',2)
plot([fzc(idx_fres) fzc(idx_fres)],[S_ms_pha(idx_fres) S_ms_pha(idx_fres)],'ko','markersize',6,'markerfacecolor','k','linew',2)

grid on
%%

figure(2),clf
subplot(121)
hold on
plot(real(S_ms),imag(S_ms),'b.')
grid on, hold on, axis square
set(gca,'xlim',[0 1],'ylim',[-.3 .7])
plot(get(gca,'xlim'),[0 0],'k--','linew',1)
plot([0 0],get(gca,'ylim'),'k--','linew',1)

plot(real(S_ms(idx_fres)),imag(S_ms(idx_fres)),'rp','markersize',4,'markerfacecolor','r','linew',2)
plot([real(S_ms(idx_fres-10)) real(S_ms(idx_fres+10))],[imag(S_ms(idx_fres-10)) imag(S_ms(idx_fres+10))],'m-','linew',1.5)
plot(real(S_ms(idx_fres+10)),imag(S_ms(idx_fres+10)),'>','markersize',4,'markeredgecolor','m','linew',1.5)
xlabel('Re(S21) [a.u.]')
ylabel('Im(S21) [a.u.]')
title('Uncalibrated')

subplot(122)
hold on
plot(real(S_ms_cal),imag(S_ms_cal),'b.')
grid on, hold on, axis square
set(gca,'ylim',[-.9e5 .9e5],'xlim',[0 1.8e5])
% plot(get(gca,'xlim'),[0 0],'k--','linew',1)
plot([0 0],get(gca,'ylim'),'k--','linew',1)
plot(get(gca,'xlim'),[0 0],'k--','linew',1)
plot(real(S_ms_cal(idx_fres)),imag(S_ms_cal(idx_fres)),'rp','markersize',4,'markerfacecolor','r','linew',2)
plot([real(S_ms_cal(idx_fres-10)) real(S_ms_cal(idx_fres+10))],[imag(S_ms_cal(idx_fres-10)) imag(S_ms_cal(idx_fres+10))],'m-','linew',1.5)
plot(real(S_ms_cal(idx_fres+10)),imag(S_ms_cal(idx_fres+10)),'^','markersize',4,'markeredgecolor','m','linew',1.5)
xlabel('Re(S21) [Hz]')
ylabel('Im(S21) [Hz]')
title('Calibrated: fres at In-phase axis')
%%
figure(3),clf
subplot(121)
hold on
plot(real(S_ms),imag(S_ms),'b.')
grid on, hold on, axis square
set(gca,'xlim',[0 1],'ylim',[-.3 .7])
plot(get(gca,'xlim'),[0 0],'k--','linew',1)
plot([0 0],get(gca,'ylim'),'k--','linew',1)

plot(real(S_ms(idx_fres)),imag(S_ms(idx_fres)),'rp','markersize',4,'markerfacecolor','r','linew',2)
plot([real(S_ms(idx_fres-10)) real(S_ms(idx_fres+10))],[imag(S_ms(idx_fres-10)) imag(S_ms(idx_fres+10))],'m-','linew',1.5)
plot(real(S_ms(idx_fres+10)),imag(S_ms(idx_fres+10)),'>','markersize',4,'markeredgecolor','m','linew',1.5)
xlabel('Re(S21) [a.u.]')
ylabel('Im(S21) [a.u.]')
title('Uncalibrated')

subplot(122)
hold on
plot(real(S_ms_cal_v),imag(S_ms_cal_v),'b.')
grid on, hold on, axis square,axis equal
set(gca,'xlim',[-.9e5 .9e5],'ylim',[-1.8e5 0])
plot(get(gca,'xlim'),[0 0],'k--','linew',1)
plot([0 0],get(gca,'ylim'),'k--','linew',1)
plot(real(S_ms_cal_v(idx_fres)),imag(S_ms_cal_v(idx_fres)),'rp','markersize',4,'markerfacecolor','r','linew',2)
plot([real(S_ms_cal_v(idx_fres-10)) real(S_ms_cal_v(idx_fres+10))],[imag(S_ms_cal_v(idx_fres-10)) imag(S_ms_cal_v(idx_fres+10))],'m-','linew',1.5)
plot(real(S_ms_cal_v(idx_fres+10)),imag(S_ms_cal_v(idx_fres+10)),'>','markersize',4,'markeredgecolor','m','linew',1.5)
xlabel('Re(S21) [Hz]')
ylabel('Im(S21) [Hz]')
title('Calibrated: fres at Quadrature axis')

%%
figure(4),clf

subplot(121)
hold on
plot(real(S_ms),imag(S_ms),'b.')
set(gca,'xlim',[0 1],'ylim',[-.3 .7])
grid on, hold on, axis square
plot(get(gca,'xlim'),[0 0],'k--','linew',1)
plot([0 0],get(gca,'ylim'),'k--','linew',1)
plot([0 real(S_ms(idx_ftone))],[0,imag(S_ms(idx_ftone))],'color',[0.8500 0.3250 0.0980],'linew',2)
plot(real(S_ms(idx_ftone)),imag(S_ms(idx_ftone)),'kx','markersize',8,'markerfacecolor','k','linew',1.5)
plot(real(S_ms(idx_fres)),imag(S_ms(idx_fres)),'rp','markersize',4,'markerfacecolor','r','linew',1.5)
xlabel('Re(S21) [a.u.]')
ylabel('Im(S21) [a.u.]')
title('Uncalibrated')

subplot(122)
hold on
plot(real(S_ms_cal),imag(S_ms_cal),'b.')
set(gca,'xlim',[-1e4 17e4],'ylim',[-.9e5 .9e5])
grid on, hold on, axis square

plot([0 0],get(gca,'ylim'),'k--','linew',1)
plot(get(gca,'xlim'),[0 0],'k--','linew',1)
plot([0 real(S_ms_cal(idx_ftone))],[0,imag(S_ms_cal(idx_ftone))],'color',[0.8500 0.3250 0.0980],'linew',2)
plot(real(S_ms_cal(idx_ftone)),imag(S_ms_cal(idx_ftone)),'kx','markersize',8,'markerfacecolor','k','linew',1.5)
plot(real(S_ms_cal(idx_fres)),imag(S_ms_cal(idx_fres)),'rp','markersize',4,'markerfacecolor','r','linew',1.5)
plot([0 0],[0,imag(S_ms_cal(idx_ftone))],'color',[0.4940 0.1840 0.5560],'linew',3)
plot([0 real(S_ms_cal(idx_fres))],[imag(S_ms_cal(idx_ftone)),imag(S_ms_cal(idx_ftone))],':','color',[0.4660 0.6740 0.1880],'linew',2)
xlabel('Re(S21) [Hz]')
ylabel('Im(S21) [Hz]')
title('Calibrated: fres at In-phase axis')

%%
% estimated freq.offset V.S. actual freq.offset
% M = 100;
% idx_frex = idx_fres-M:idx_fres+M;
% Df_est_arr = imag(S_ms_cal(idx_frex));
Df_est_arr = imag(S_ms_cal);

figure(5),clf
hold on
% axis square
plot(fzc./BW,Df_est_arr./BW,'b-','linew',2)
plot(fzc./BW,fzc./BW,'k-','linew',2)
legend({'frequency offset estimate','true frequency offset'},'location','northwest')
grid on
xlim([-1 1])
ylim([-1 1])
xlabel('Frequency Offset [fractional Bandwidth]')
ylabel('Frequency error estimation [fractional Bandwidth]')