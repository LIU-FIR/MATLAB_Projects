%
% Demonstrate calibration,tone tracking and flux-ramp demodulation, with
% real resonator's S21 parameters measure by GuoFu(Roach2).
%
clear;clc;close all;

%%
% S21 parameters, Qr = 1/(1/Qi+1/Qc)
% 
fres = 4.68e9; % original resonance frequency
Qr = 5858.96;
Qc = 23480.66;
phase_dly = 0.387; % system phase delay.
% computing parameters
BW = 100e3; % Yu SLAC Microresonator RF (SMuRF) Electronics, P8, CMB umux resonators bandwidth
fstep = 5e1; % computational frequency interval/step.
f = fres-BW/2*10:fstep:fres+BW/2*10; % computational frequency range, fres is central frequency. 
fzc = f-fres; % move fres to zero frequency,
%
S_ms = S21_Rp(Qr,Qc,fres,f-fres,phase_dly); % measured S21 with phase delay.
S_ms_pha = angle(S_ms);
idx_fres = find(abs(S_ms)==min(abs(S_ms)));
%%
% S21模型
figure(1),clf
hold on
yyaxis left
plot(fzc,abs(S_ms),'b-','linew',2)
%
plot([fzc(idx_fres) fzc(idx_fres)],[abs(S_ms(idx_fres)) abs(S_ms(idx_fres))],'ko','markersize',6,'markerfacecolor','k','linew',2)
xlabel('Hz')
ylabel('Magnitude: (a.u.)')
%
yyaxis right
plot(fzc,S_ms_pha,'r-','linew',2)
ylabel('Phase: (rad)')
yrange = get(gca,'ylim');
xrange = get(gca,'xlim');
plot([fzc(idx_fres) fzc(idx_fres)],[yrange(1) yrange(2)],'k--','linew',2)
plot([fzc(idx_fres) fzc(idx_fres)],[S_ms_pha(idx_fres) S_ms_pha(idx_fres)],'ko','markersize',6,'markerfacecolor','k','linew',2)
%%
Eta = Eta_gen(S_ms,fres,f,'real');
disp(['Eta = ',num2str(Eta)]);
S_ms_calib = Eta.*S_ms;
%%
figure(2),clf
subplot(121)
hold on
plot(real(S_ms),imag(S_ms),'b.')
grid on, hold on, axis square
set(gca,'xlim',[0.4 1],'ylim',[-.1 .5])
plot(get(gca,'xlim'),[0 0],'k--','linew',1)
plot([0 0],get(gca,'ylim'),'k--','linew',1)

plot(real(S_ms(idx_fres)),imag(S_ms(idx_fres)),'rp','markersize',4,'markerfacecolor','r','linew',2)
% plot([real(S_ms(idx_fres-10)) real(S_ms(idx_fres+10))],[imag(S_ms(idx_fres-10)) imag(S_ms(idx_fres+10))],'m-','linew',1.5)
% plot(real(S_ms(idx_fres+10)),imag(S_ms(idx_fres+10)),'>','markersize',4,'markeredgecolor','m','linew',1.5)
xlabel('Re(S21) [a.u.]')
ylabel('Im(S21) [a.u.]')
title('Uncalibrated')

subplot(122)
hold on
plot(real(S_ms_calib),imag(S_ms_calib),'b.')
grid on, hold on, axis square
set(gca,'ylim',[-3e5 3e5],'xlim',[9e5 15e5])
% plot(get(gca,'xlim'),[0 0],'k--','linew',1)
plot([0 0],get(gca,'ylim'),'k--','linew',1)
plot(get(gca,'xlim'),[0 0],'k--','linew',1)
plot(real(S_ms_calib(idx_fres)),imag(S_ms_calib(idx_fres)),'rp','markersize',4,'markerfacecolor','r','linew',2)
% plot([real(S_ms_cal(idx_fres-10)) real(S_ms_cal(idx_fres+10))],[imag(S_ms_cal(idx_fres-10)) imag(S_ms_cal(idx_fres+10))],'m-','linew',1.5)
% plot(real(S_ms_cal(idx_fres+10)),imag(S_ms_cal(idx_fres+10)),'^','markersize',4,'markeredgecolor','m','linew',1.5)
xlabel('Re(S21) [Hz]')
ylabel('Im(S21) [Hz]')
title('Calibrated: fres at In-phase axis')

%%
n_quanta = 4; % flux quanta number.
framp = 10e3; % sawtooth rate.
srate = 1e6; % sampling flux ramp rate
fc  = framp * n_quanta; % the fundemental frequency determined by flux ramp modulation.
n_horder = 3; % use n_harmonics to approximate modulated fres.
delay  = 0; % defined by system delay (filter, cryo, etc.) 
% mu = 5e-2; % user defined closed-loop gain
mu = 8e-2;
%
lambda = 0.3;
n      = 0:2^12-1; % all recorded samples
tn      = n/srate; % sample interval (secs)
fpp = 100e3; % squid frequency swing, Yu,2022, P5

nT_det1 = (srate/framp)*2; % 接收到光子1的时刻索引
ndT_det1 = (srate/framp)*1; % 接收到光子1持续时间索引
ph1 = pi/20; % detector signal inducing phase change1.

nT_det2 = (srate/framp)*3; % 接收到光子2的时刻索引
ndT_det2 = (srate/framp)*1; % 接收到光子2持续时间索引
ph2 = pi/15; % detector signal inducing phase change1.

nT_det3 = (srate/framp)*4; % 接收到光子2的时刻索引
ndT_det3 = (srate/framp)*1; % 接收到光子2持续时间索引
ph3 = pi/30; % detector signal inducing phase change1.


D_squid_free  = lambda*sin(2*pi*fc*tn)./(1 + lambda*sin(2*pi*fc*tn));
D_squid_1  = lambda*sin(2*pi*fc*tn + ph1)./(1 + lambda*sin(2*pi*fc*tn + ph1));
D_squid_2 = lambda*sin(2*pi*fc*tn + ph2)./(1 + lambda*sin(2*pi*fc*tn + ph2));
D_squid_3 = lambda*sin(2*pi*fc*tn + ph3)./(1 + lambda*sin(2*pi*fc*tn + ph3));


% D_squid_free  = lambda*cos(2*pi*fc*tn)./(1 + lambda*cos(2*pi*fc*tn));
% D_squid_1  = lambda*cos(2*pi*fc*tn + ph1)./(1 + lambda*cos(2*pi*fc*tn + ph1));
% D_squid_2 = lambda*cos(2*pi*fc*tn + ph2)./(1 + lambda*cos(2*pi*fc*tn + ph2));
% D_squid_3 = lambda*cos(2*pi*fc*tn + ph3)./(1 + lambda*cos(2*pi*fc*tn + ph3));

phases = @(a,n) a*ones(1,numel(n));
phases_det = phases(0,n);

phases_det(nT_det1+1:nT_det1+ndT_det1) = phases(ph1,nT_det1+1:nT_det1+ndT_det1);
phases_det(nT_det1+ndT_det1+1:nT_det1+2*ndT_det1) = phases(ph2,nT_det1+ndT_det1+1:nT_det1+2*ndT_det1);
phases_det(nT_det1+2*ndT_det1+1:nT_det1+3*ndT_det1) = phases(ph3,nT_det1+2*ndT_det1+1:nT_det1+3*ndT_det1);
            
D_squid_det = [D_squid_free(1:nT_det1),...
                D_squid_1(nT_det1+1:nT_det1+ndT_det1),...
                D_squid_2(nT_det1+ndT_det1+1:nT_det1+2*ndT_det1),...
                D_squid_3(nT_det1+2*ndT_det1+1:nT_det1+3*ndT_det1),...
                D_squid_free(nT_det1+3*ndT_det1+1:end)];            
            
B = fpp./(max(D_squid_free)-min(D_squid_free));
f_squid_free = B*D_squid_free;
f_squid = B*D_squid_det;
flux_ramp = B/3*sawtooth(2*pi*framp*tn);

%%
figure(3),clf
subplot(211)
hold on
yyaxis left
plot(f_squid,'r-','linew',2)
plot(f_squid_free,'k-','linew',2)
xlabel('sample pts')
ylabel('Squid frequency (Hz)')
% plot(flux_ramp,'m-','linew',2)
yyaxis right
plot(phases_det,'b-','linew',3)
xlim([0 1500])
xlabel('sample pts')
ylabel('Phases (rad)')
subplot(212)
hold on
plot(flux_ramp,'m-','linew',2)
xlim([0 1500])
%%
% LMS filter
M = 2*n_horder+1;
alpha = zeros(M,1);
fy = zeros(1,numel(n)); % LMS filter output
theta_cal = zeros(1,numel(n));

theta_cal_n = zeros(1,numel(n));
est_foffset = zeros(1,numel(n));

for i = 0:numel(n)-1
    idx = i+1; % array index
    s_i = harmonics_gen(n_horder,fc,srate,i);% current time's harmonics.
    fy(idx) = dot(alpha,s_i);  % LMS filter's output  
    
    % Compute demodulated phase(dector)
    % when f_squid = lambda*sin(2*pi*fc*tn)./(1 + lambda*sin(2*pi*fc*tn)); 
    theta_cal(idx) = atan2(alpha(2),alpha(1)); 
    % when f_squid = lambda*cos(2*pi*fc*tn)./(1 + lambda*cos(2*pi*fc*tn)); 
    % theta_cal(idx) = atan2(-alpha(1),alpha(2)); 

    e_i_est = imag(S21_Rp(Qr,Qc,fres,f_squid(idx) - fy(idx),phase_dly)*Eta);
    
    alpha = alpha + mu*e_i_est * s_i'; % update filter coefficients. 
    est_foffset(idx) = imag(S21_Rp(Qr,Qc,fres,f_squid(idx),phase_dly)*Eta);
end
%%
figure(4),clf
subplot(211)
hold on
% plot(tn,f_squid,'k-','linew',1)
plot(f_squid,'ks-','linew',1)
plot(est_foffset,'m.','linew',1)
% plot(tn,fy,'b.-','linew',1)
% plot(fy,'b.-','linew',1)
% plot([T_det,T_det],[-max(f_squid),max(f_squid)],'r-','linew',1)

xlabel('Sampling pts')
ylabel('Tracking SQUID frequency [Hz]')
% xlim([0 1500]./srate)
xlim([0 100])

subplot(212)
% plot(tn,f_squid-fy,'r-','linew',1)
plot(f_squid-fy,'r-','linew',1)
xlabel('Sampling pts')
ylabel('Frequency error [Hz]')
% xlim([0 1500]./srate)
xlim([0 1500])

%
err_relative = (theta_cal-phases_det)./phases_det;

%%
figure(5),clf
subplot(311)
hold on
plot(f_squid,'k-','linew',2)
plot(fy,'b.-','linew',1)
xlim([0 1500])
xlabel('sample pts')
ylabel('Squid frequency (Hz)')
subplot(312)
hold on
plot(flux_ramp,'m-','linew',2)
xlim([0 1500])
xlabel('sample pts')
ylabel('Flux ramp sawtooth amplitude (a.u.)')
subplot(313)
hold on
plot(phases_det,'k-','linew',2)
plot(theta_cal,'b-','linew',2)
ylim([-.1 .3])
xlim([0 1500])
xlabel('sample pts')
ylabel('Phases (Hz)')

%%
% 利用FFT后的频域复相关(Hadmard乘积)计算相位
n1 = n(1:100);
f_squid_free = est_foffset(1:100);
f_squid_ph1 = est_foffset(nT_det1+1:nT_det1+ndT_det1);
f_squid_ph2 = est_foffset(nT_det1+ndT_det1+1:nT_det1+2*ndT_det1);
f_squid_ph3 = est_foffset(nT_det1+2*ndT_det1+1:nT_det1+3*ndT_det1);
F_squid_free = fft(f_squid_free);
F_squid_ph1 = fft(f_squid_ph1);
F_squid_ph2 = fft(f_squid_ph2);
F_squid_ph3 = fft(f_squid_ph3);

XC1 = F_squid_ph1 .* conj(F_squid_free);
XC2 = F_squid_ph2 .* conj(F_squid_free);
XC3 = F_squid_ph3 .* conj(F_squid_free);

N = length(n1);
% ,'markerfacecol','k'
figure(7),clf
subplot(231)
hold on
yyaxis left
% stem(0:N-1,abs(X),'ko-','markersize',2)
stem(0:N-1,abs(XC1),'bs-','markerfacecol','b','markersize',2)
xlim([0,N/4])
ylabel('Magnitude (a.u)')
yyaxis right
plot(0:N-1,angle(XC1),'ro--','markerfacecol','r','markersize',2)
xlim([0,N/4])
xlabel('FFT bins')
ylabel('Phases (rad)')

subplot(232)
hold on
yyaxis left
% stem(0:N-1,abs(X),'ko-','markersize',2)
stem(0:N-1,abs(XC2),'bs-','markerfacecol','b','markersize',2)
xlim([0,N/4])
ylabel('Magnitude (a.u)')
yyaxis right
plot(0:N-1,angle(XC2),'ro--','markerfacecol','r','markersize',2)
xlim([0,N/4])
xlabel('FFT bins')
ylabel('Phases (rad)')

subplot(233)
hold on
yyaxis left
% stem(0:N-1,abs(X),'ko-','markersize',2)
stem(0:N-1,abs(XC3),'bs-','markerfacecol','b','markersize',2)
xlim([0,N/4])
ylabel('Magnitude (a.u)')
yyaxis right
plot(0:N-1,angle(XC3),'ro--','markerfacecol','r','markersize',2)
xlim([0,N/4])
xlabel('FFT bins')
ylabel('Phases (rad)')

subplot(212)
hold on
plot(phases_det,'k-','linew',2)
plot(nT_det1+1:nT_det1+ndT_det1,angle(XC1(5))*ones(1,100),'rs','markerfacecol','r','markersize',1.2)
plot(nT_det1+ndT_det1+1:nT_det1+2*ndT_det1,angle(XC2(5))*ones(1,100),'rs','markerfacecol','r','markersize',1.2)
plot(nT_det1+2*ndT_det1+1:nT_det1+3*ndT_det1,angle(XC3(5))*ones(1,100),'rs','markerfacecol','r','markersize',1.2)
xlim([0 800])
ylim([0 0.25])
xlabel('Samples (pts)')
ylabel('Phases (rad)')
%
% compute the relative error
disp([(angle(XC1(5))-ph1)/ph1, (angle(XC2(5))-ph2)/ph2, (angle(XC3(5))-ph3)/ph3])

%%
% 利用时域拟合求相位移动，注意用频率offset的估计值
n1 = n(1:100);
f_squid_free = est_foffset(1:100);
f_squid_ph1 = est_foffset(nT_det1+1:nT_det1+ndT_det1);
f_squid_ph2 = est_foffset(nT_det1+ndT_det1+1:nT_det1+2*ndT_det1);
f_squid_ph3 = est_foffset(nT_det1+2*ndT_det1+1:nT_det1+3*ndT_det1);

fit_free = fit(n1',f_squid_free','fourier3');
fit_ph1 = fit(n1',f_squid_ph1','fourier3');
fit_ph2 = fit(n1',f_squid_ph2','fourier3');
fit_ph3 = fit(n1',f_squid_ph3','fourier3');

% T1_fit = 2*pi/fit_free.w;
% fc1_fit = srate/(T1_fit);
% wc1_fit = 2*pi*fc1_fit/srate;

w_free = fit_free.w;
w_ph1 = fit_ph1.w;
w_ph2 = fit_ph2.w;
w_ph3 = fit_ph3.w;

f_squid_free_fit = fit_free.a0 + fit_free.a1*cos(w_free*n1) + fit_free.b1*sin(w_free*n1)...
    + fit_free.a2*cos(2*w_free*n1) + fit_free.b2*sin(2*w_free*n1) ...
    + fit_free.a3*cos(3*w_free*n1) + fit_free.b3*sin(3*w_free*n1);

f_squid_ph1_fit = fit_ph1.a0 + fit_ph1.a1*cos(w_ph1*n1) + fit_ph1.b1*sin(w_ph1*n1)...
    + fit_ph1.a2*cos(2*w_ph1*n1) + fit_ph1.b2*sin(2*w_ph1*n1) ...
    + fit_ph1.a3*cos(3*w_ph1*n1) + fit_ph1.b3*sin(3*w_ph1*n1);

f_squid_ph2_fit = fit_ph2.a0 + fit_ph2.a1*cos(w_ph2*n1) + fit_ph2.b1*sin(w_ph2*n1)...
    + fit_ph2.a2*cos(2*w_ph2*n1) + fit_ph2.b2*sin(2*w_ph2*n1) ...
    + fit_ph2.a3*cos(3*w_ph2*n1) + fit_ph2.b3*sin(3*w_ph2*n1);

f_squid_ph3_fit = fit_ph3.a0 + fit_ph3.a1*cos(w_ph3*n1) + fit_ph3.b1*sin(w_ph3*n1)...
    + fit_ph3.a2*cos(2*w_ph3*n1) + fit_ph3.b2*sin(2*w_ph3*n1) ...
    + fit_ph3.a3*cos(3*w_ph3*n1) + fit_ph3.b3*sin(3*w_ph3*n1);

phi1 = atan2(fit_ph1.b1,fit_ph1.a1) - atan2(fit_free.b1,fit_free.a1);
phi2 = atan2(fit_ph2.b1,fit_ph2.a1) - atan2(fit_free.b1,fit_free.a1);
phi3 = atan2(fit_ph3.b1,fit_ph3.a1) - atan2(fit_free.b1,fit_free.a1);

% compute the relative error
disp([(-phi1-ph1)/ph1, (-phi2-ph2)/ph2, (-phi3-ph3)/ph3])
%%
figure(6),clf
subplot(231)
hold on
plot(n1,f_squid_free,'ko','markerfacecol','k','markersize',1.5)
plot(n1,f_squid_free_fit,'k-')
plot(n1,f_squid_ph1,'rs','markerfacecol','r','markersize',1.5)
plot(n1,f_squid_ph1_fit,'r-')
xlim([0 25])
ylabel('Squid frequency (Hz)')
xlabel('Samples (pts)')
subplot(232)
hold on
plot(n1,f_squid_free,'ko','markerfacecol','k','markersize',1.5)
plot(n1,f_squid_free_fit,'k-')
plot(n1,f_squid_ph2,'rs','markerfacecol','r','markersize',1.5)
plot(n1,f_squid_ph2_fit,'r-')
xlim([0 25])
ylabel('Squid frequency (Hz)')
xlabel('Samples (pts)')
subplot(233)
hold on
plot(n1,f_squid_free,'ko','markerfacecol','k','markersize',1.5)
plot(n1,f_squid_free_fit,'k-')
plot(n1,f_squid_ph3,'rs','markerfacecol','r','markersize',1.5)
plot(n1,f_squid_ph3_fit,'r-')
xlim([0 25])
ylabel('Squid frequency (Hz)')
xlabel('Samples (pts)')
subplot(212)
hold on
plot(phases_det,'k-','linew',2)
plot(nT_det1+1:nT_det1+ndT_det1,-phi1*ones(1,100),'rs','markerfacecol','r','markersize',1.2)
plot(nT_det1+ndT_det1+1:nT_det1+2*ndT_det1,-phi2*ones(1,100),'rs','markerfacecol','r','markersize',1.2)
plot(nT_det1+2*ndT_det1+1:nT_det1+3*ndT_det1,-phi3*ones(1,100),'rs','markerfacecol','r','markersize',1.2)
xlim([0 800])
ylim([0 0.25])
xlabel('Samples (pts)')
ylabel('Phases (rad)')
%%
