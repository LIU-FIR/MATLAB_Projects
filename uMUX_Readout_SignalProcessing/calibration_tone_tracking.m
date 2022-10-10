clear;clc;close all;

%%
Qi = 348000;
df1 = 10e3; % Hz % asymmetry, resonance frequency deviation.
df = 0;
% Qr = 1/(1/Qi+1/Qc); % relationship between Qr and Qi,Qc.
Qr = 17000;
Qc = 20000;
fres = 5e9; % original resonance frequency


BW = 100e3; % Yu SLAC Microresonator RF (SMuRF) Electronics, P8, CMB umux resonators bandwidth
fstep = 5e1; % computational frequency interval/step.
f = fres-BW/2*10:fstep:fres+BW/2*10; % computational frequency range, fres is central frequency. 
fzc = f-fres; % move fres to zero frequency,
%
S =1-Qr/Qc*1./(1+1i*2*Qr*(f-fres)./fres); % Ideal S21, fres at real axis.
% S_ms = 1-Qr/Qc*(1-1i*2*Qc*df1/fr)./(1+1i*2*Qr*(f-fr)./fr); % S21 measurement,asymmetry.
S_ms = S.*exp(1i*pi/8); % uncalibrated S21 measurement，fres is not located at real axis.

idx_offset = 10;
idx_fres = find(abs(S_ms)==min(abs(S_ms))); % "shifted" resonance frequency.
idx_ftone = idx_fres + idx_offset; %% assume ftone stays unchanged.

%
Eta = (fzc(idx_fres+1)-fzc(idx_fres-1))./(S_ms(idx_fres+1)-S_ms(idx_fres-1));
% Eta1 = (fk(idx_fres+1)-fk(idx_fres-1))./(S_ms(idx_fres-1)-S_ms(idx_fres+1));
disp(['Eta = ',num2str(Eta)]);
S_ms_cal = Eta.*S_ms*1i;% calibrated S21 to vertical axis, Eta is a rotated vector.
S_ms_cal_v = Eta.*S_ms; % calibrated S21 to horizontal axis, Eta is a rotated vector.

Df_est = imag(S_ms_cal(idx_ftone));
Df_est1 = imag((1-Qr/Qc*1./(1+1i*2*Qr*(fzc(idx_ftone))./fres)).*exp(1i*pi/8)*Eta*1i);
disp([Df_est fzc(idx_ftone) Df_est1])

fres_offset = -fstep*100;
S_ms_fresOffset = (1-Qr/Qc*1./(1+1i*2*Qr*((f-fres_offset)-fres)./fres)).*exp(1i*pi/8);
figure(1),clf
hold on
% plot(fzc,abs(S_ms),'b.-','linew',1)
plot(fzc,abs(S_ms_fresOffset),'m.-','linew',1)
plot(fzc(find(f==fres)),abs(S_ms_fresOffset(find(f==fres))),'kx')
xlim([-2e5 2e5])

%%
srate = 2.4e6; % sampling rate
flms  = 16e3; % the fundemental frequency determined by flux ramp modulation.
n_horder = 3; % use n_harmonics to approximate modulated fres.
delay  = 0; % defined by system delay (filter, cryo, etc.) 
mu = 1e-2; % user defined closed-loop gain
%
lambda = 0.3;
n      = 0:2^14-1;
tn      = n/srate;
fpp = 100e3; % Yu,2022, P5
% D_squid  = lambda*sin(2*pi*(flms/srate)*n)./(1 + lambda*sin(2*pi*(flms/srate)*n));
D_squid  = lambda*sin(2*pi*flms*tn)./(1 + lambda*sin(2*pi*flms*tn));
B = fpp./(max(D_squid)-min(D_squid));
f_squid = B*D_squid;

% figure(2),clf
% plot(f_squid,'k-','linew',1)
% xlim([0 1500])
%%
% LMS filter
M = 2*n_horder+1;
alpha = zeros(M,1);
fy = zeros(1,numel(n)); % LMS filter output


for i = 0:numel(n)-1
    idx = i+1; % array index
    s_i = harmonics_gen(n_horder,flms,srate,i);
%     if idx==1
%         y(idx)=imag(S_ms_cal(idx_ftone));
%     else
%         y(idx) = dot(alpha,s_i);
%     end
    fy(idx) = dot(alpha,s_i);    
%     idx_f = dsearchn(fzc',y(idx));
    e_i = imag((1-Qr/Qc*1./(1+1i*2*Qr*(f_squid(idx))./fres)).*exp(1i*pi/8)*Eta*1i)- fy(idx);
%     e_i = f_squid(idx) - fy(idx);
    alpha = alpha + mu*e_i*s_i';   
end
%%
err = f_squid-fy;
err_rms = rms(err(2000:end));
err_mean = mean(err(2000:end));
figure(2),clf
subplot(211)
hold on
plot(f_squid,'k-','linew',1)
plot(fy,'b.-','linew',1)
xlabel('Sampling pts')
ylabel('Tracking SQUID frequency [Hz]')
xlim([0 2000])
subplot(212)
plot(f_squid-fy,'r-','linew',1)
xlabel('Sampling pts')
ylabel('Frequency error [Hz]')
xlim([0 2000])


%%
figure(3),clf
plot(f_squid,'k-','linew',1)
xlabel('Sampling pts')
ylabel('SQUID frequency [Hz]')
xlim([0 1500])
