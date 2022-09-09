%%
% flux ramp modulation的模拟
% 参考文献：The The Microwave SQUID Multiplexer_PHD_Thesis(Mates) 
%%
clear; clc;
phi_0 = 1e-15; % 超导量子磁通，典型值为 10^(-15)Wb
lamb = 1/3; % The The Microwave SQUID Multiplexer_PHD_Thesis, Ch2.1.2
Lc = 77.6e-12;% 77.6 pH, Ch2.1.3, fig.2.6
Mc = 1.65e-12;% 1.65 pH, Ch2.1.3, fig.2.6
Ls = 18.9e-12;% 18.9 pH, Ch2.1.3, fig.2.6

phi_wave = linspace(-phi_0,phi_0,200);

L_phi = Lc-Mc^2/Ls*(lamb*cos(2*pi*phi_wave/phi_0)./(1+lamb*cos(2*pi*phi_wave/phi_0)));

%%
figure(1)
plot(phi_wave,L_phi/1e-12,'k-')
ylim([77.5 77.7])
grid on
%%
clear;clc;
f_ramp= 5e3; % sawtooth ramp频率, Ch6.3.4
n_phi0 = 10; % Ch6.3.4
phi_i = pi/10;
f_RF = 1.024e9;% SQUID的RF信号
% f_RF = 1.024e6;% SQUID的RF信号
fs = 4.096e9;
% fs = 4.096e6;
N = 40000;
n = 0:N-1;

theta_i = cos(2*pi*n_phi0*f_ramp/fs*n + phi_i);
r_bias = 0.2+0.3i; % 引入信号的偏置项
r = exp(1i*(2*pi*f_RF/fs*n + theta_i))+r_bias; % SQUID的响应



%两种方法反求调制相位项的结果相同
LO1 = exp(-1i*2*pi*f_RF/fs*n);
r_mix = r.*LO1;
% theta_i_cal = atan2(imag(r_mix),real(r_mix)); % 通过混频消除f_RF引起的快速相位变化

theta_i_cal = unwrap(atan2(imag(r),real(r)))-2*pi*f_RF/fs*n;% 直接通过I,Q来求展开后的相位，然后扣除每一时刻的f_RF引起相位变化
theta_i_cal_unbiased = unwrap(atan2(imag(r-r_bias),real(r-r_bias)))-2*pi*f_RF/fs*n;% 直接通过I,Q来求展开后的相位，然后扣除每一时刻的f_RF引起相位变化

%%
figure(2)
subplot(2,1,1)
plot(n,theta_i,'k-','LineWidth',.5)
hold on
plot(n,theta_i_cal_unbiased,'b.')
xlim([0 50000])
subplot(2,1,2)
plot(n,theta_i_cal,'r-')
xlim([0 50000])
%%
figure(3) % 画相图
% scatter(real(r_mix),imag(r_mix))
scatter(real(r)-0.2,imag(r)-0.3,'b.')
hold on
scatter(real(r-r_bias)-0.2,imag(r-r_bias)-0.3,'r.')
axis equal
%%
N1 = 100; N2 = 1000; N3 = 10000;
n1 = 0:N1-1;
n2 = 0:N2-1;
n3 = 0:N3-1;
n4 = 0:N-1;

w_RF = 2*pi*n_phi0*f_ramp/fs;
phi_i1 = atan2(-sum(theta_i_cal_unbiased(n1+1).*sin(w_RF*n1)),sum(theta_i_cal_unbiased(n1+1).*cos(w_RF*n1)));
phi_i2 = atan2(-sum(theta_i_cal_unbiased(n2+1).*sin(w_RF*n2)),sum(theta_i_cal_unbiased(n2+1).*cos(w_RF*n2)));
phi_i3 = atan2(-sum(theta_i_cal_unbiased(n3+1).*sin(w_RF*n3)),sum(theta_i_cal_unbiased(n3+1).*cos(w_RF*n3)));
phi_i4 = atan2(-sum(theta_i_cal_unbiased(n4+1).*sin(w_RF*n4)),sum(theta_i_cal_unbiased(n4+1).*cos(w_RF*n4)));

%% 2021/10/27
%参数在ShenMengping的Excel文档基础上设定

clc;clear;
Nx=2^10;Ny=2^8;
N = Nx*Ny; % Buffer长度
n = 0:N-1;
fs=4.096e9; %采样率


f_RF = 2*fs/Nx; % 射频信号频率
f_ramp = 500e3; % sawtooth ramp频率, Becker 2019
n_phi0 = 2; % Becker 2019
fc = n_phi0*f_ramp; % 调制相位的变化频率

A = 1;
B = 1;

w_RF = 2*pi*f_RF/fs;
wc = 2*pi*fc/fs;

phi_i = pi/30; % 光子信号通过TES电路感应出来的相位
theta_i = B*cos(wc*n+phi_i); % 经过flux ramp调制的相位

% x = A*cos(w_RF*n + B*sin(wc*n));
s_IQ = A*exp(1i*(w_RF*n + theta_i));
SNR_dB = 30; % 信噪比 (dB)
SNR_mag = db2mag(SNR_dB);

Aeff = 2^0.5*A/2; % 正弦信号的有效值
noise1 = Aeff/SNR_mag * complex(randn(1,N),randn(1,N));
% noise1 = 0;
x_IQ = s_IQ + noise1;
x_IQ_real = real(x_IQ);

theta_i_x = unwrap(atan2(imag(x_IQ),real(x_IQ))); %根据实际信号x(带噪声)计算的总相位
theta_i_x_mod = unwrap(atan2(imag(x_IQ),real(x_IQ))) - w_RF*n; %根据实际信号x(带噪声)计算的调制相位
%%
figure(1)
plot(n,real(x_IQ),'k-')
xlim([0 60000])
ylim([-1.5*A 1.5*A])
grid on

figure(2)
subplot(2,1,1)
plot(n,theta_i_x,'k-','LineWidth',1)
subplot(2,1,2)
plot(n,theta_i_x_mod,'r-','LineWidth',1)
hold on
plot(n,theta_i,'b-','LineWidth',2)
xlim([0 length(n)/10])
grid on
%%
Nspr = fs/fc;
n1 = 0:Nspr-1;
phi_i_cal = atan(-sum(theta_i(1:Nspr).*sin(wc*n1))/sum(theta_i(1:Nspr).*cos(wc*n1)));
phi_i_x_cal = atan(-sum(theta_i_x_mod(1:Nspr).*sin(wc*n1))/sum(theta_i_x_mod(1:Nspr).*cos(wc*n1)));
r_error = abs(phi_i_x_cal - phi_i)/phi_i;
display(r_error)
%% 2022-05-16
clc;clear;
Nx=2^10;Ny=2^9;
N = Nx*Ny; % DAC Buffer-length
n = 0:N-1;
srate_DAC = 4.096e9; % DAC sampling rate
srate_ADC = 512e6; % ADC sampling rate
Rd = srate_DAC/srate_ADC; 

f1 = 8e6;
f2 = 128e6;
f3 = 144e6;
f4 = 160e6;
fc = 15.625e3; % modulate phase's frequency
A = 8190/4;

phi_wave = sin(2*pi*fc/srate_DAC*n);
x = A*cos(2*pi*f1/srate_DAC*n + phi_wave)...
    + A*cos(2*pi*f2/srate_DAC*n)...
    + A*cos(2*pi*f3/srate_DAC*n + phi_wave)...
    + A*cos(2*pi*f4/srate_DAC*n);
%%
figure(1),clf
subplot(211)
plot(n,x,'b-')
xlim([0 3200])
%%
y = downsample(x,Rd);
m = 0:N/Rd-1;

subplot(212)
plot(m,y,'k-')
xlim([0 400])
%%
Nfft = 512;
y_win = y(1:Nfft);
Y = abs(fft(y_win));
k = (-Nfft/2:Nfft/2-1);
fk = (-Nfft/2:Nfft/2-1)*srate_ADC/Nfft/1e6; % MHz

figure(2),clf
% plot(fk,fftshift(Y),'k-')
plot(k,fftshift(Y),'k-')
%%
Ly = length(y); 
N_frames = Ly/Nfft;
srate_fft = srate_ADC/Nfft; % the sample rate at FFT's output.

y_frames = reshape(y,Nfft,N_frames);
Y_frames = fft(y_frames);

idx_frex = [8+1 128+1 144+1 160+1];
Y_f1 = Y_frames(idx_frex(1),:); % 8MHz-bin, phase-modulated
Y_f2 = Y_frames(idx_frex(2),:);
Y_f3 = Y_frames(idx_frex(3),:); % 144MHz-bin, phase-modulated
Y_f4 = Y_frames(idx_frex(4),:);

figure(3),clf
hold on
plot(angle(Y_f1),'ko-')
plot(angle(Y_f3)+.1,'bs-')
plot(angle(Y_f2),'r-')
plot(angle(Y_f4)+.1,'m-')

%%
% clear;clc;
filename = 'c:\Jianguoyun\WorkingProjects\TES读出系统\MATLAB\flux_ramp_mod\aftertomb_FIXCLK_drawdata.txt';
fd = fopen(filename);
formatSpec = '%d';

Ns = 512*1024*50;
C_text = textscan(fd,formatSpec,'Delimiter','');
data = C_text{1,1};
fclose(fd);
%%
figure(4)
plot(data(1:1000),'k-')
xlim([0 400])
%%
data_win = data(1:65536*8);
data_frames = reshape(data_win,Nfft,[]);
Z_frames = fft(data_frames);

idx_frex = [8+1 128+1 144+1 160+1];
Z_f1 = Z_frames(idx_frex(1),:); % 8MHz-bin, phase-modulated
Z_f2 = Z_frames(idx_frex(2),:);
Z_f3 = Z_frames(idx_frex(3),:); % 144MHz-bin, phase-modulated
Z_f4 = Z_frames(idx_frex(4),:);
%%
figure(5)
hold on
plot(angle(Z_f1),'k-')
plot(angle(Z_f3),'b-')
plot(angle(Z_f2),'r-')
plot(angle(Z_f4),'m-')
%%
figure(6)
plot(k,fftshift(abs(Z_frames(:,1))),'b-')
%%
angle_f1 = angle(Z_f1);
angle_f1_len = length(angle_f1);

A_f1 = fft(angle_f1);
% srate_angle_fft = srate_ADC/Nfft;
ka = (-angle_f1_len/2:angle_f1_len/2-1)*srate_fft/angle_f1_len;
figure(7)
plot(ka,fftshift(abs(A_f1)),'k-')
xlim([0 500000])