%% 测试沈梦萍ADC数据文件

%% 2021/11/8
% 模拟DAC产生的数据buffer

clc;clear;
N=2^10; % FFT长度为1024
M=2^10; % 帧数
L = N*M; % Buffer长度
n = 0:L-1;
fs = 4.096e9; %采样率
dfs = 0; % 采样率偏移
fs = fs + dfs;

ks_RF = [2,32,36,40];
dk = 0.012; %  射频信号频率对应fs的FFT-bin
ks_RF_bias = ks_RF+dk;
f_RF = ks_RF_bias*fs/N; % 射频信号频率
% f_ramp = 500e3; % sawtooth ramp频率, Becker 2019
% n_phi0 = 2; % Becker 2019
% fc = 6.25e3; % 调制相位的变化频率
fc = 15625;% 调制相位的变化频率

A = 8190/4;
B = 1;

w_RF = 2*pi*f_RF/fs;
wc = 2*pi*fc/fs;

% phi_i = 0; % 光子信号通过TES电路感应出来的相位
% theta_i = B*sin(wc*n+phi_i); % 经过flux ramp调制的相位
theta_i = -abs(sin(wc*n)); % 经过flux ramp调制的相位

s1 = A*cos(w_RF(1)*n + theta_i);
s2 = A*cos(w_RF(2)*n + theta_i);
s3 = A*cos(w_RF(3)*n +theta_i);
s4 = A*cos(w_RF(4)*n + theta_i);
si = s1 + s2 + s3 + s4;
% sq = A*sin(w_RF*n + theta_i);
% s_IQ = A*exp(1i*(w_RF*n + theta_i));

si_arr = reshape(si,[N,M]);
y_arr = fft(si_arr);
yk_RF_seq = y_arr(ks_RF+1,:); 
% yk_RF_mag = abs(yk_RF_seq);
yk_RF_phase = atan2(imag(yk_RF_seq),real(yk_RF_seq));
% theta_i_x = unwrap(atan2(imag(x_IQ),real(x_IQ))); %根据实际信号x(带噪声)计算的总相位
% theta_i_x_mod = unwrap(atan2(imag(x_IQ),real(x_IQ))) - w_RF*n; %根据实际信号x(带噪声)计算的调制相位
%%
kax = -N/2:N/2-1;
figure
stem(kax,fftshift(abs(y_arr(:,1))),'k.-')
xlim([0 50])
grid on
%%
slope = 2*pi*dk;
slope_fit = polyfit(0:M-1,unwrap(yk_RF_phase(1,:)),1);
(slope - slope_fit(1))/slope
%%
theta_dk = slope*(0:M-1);
theta_dk_fit = slope_fit(1)*(0:M-1);
figure
plot(unwrap(yk_RF_phase(1,:)-theta_dk),'r-')
hold on
plot(unwrap(yk_RF_phase(1,:)),'k-')
hold on
plot(detrend(unwrap(yk_RF_phase(1,:))-theta_dk_fit),'b-')
grid on
ylim([-2 2])
%%
% 确定相位的变化频率
Np = length(yk_RF_phase);
kpha =  -M/2:M/2-1;

Zk_phase = fft(yk_RF_phase,[],2); % 沿行FFT，得到相位的变化
figure
subplot(2,2,1)
stem(-M/2:M/2-1,fftshift(abs(Zk_phase(1,:))),'k.');
xlim([-20 20])
grid on
subplot(2,2,2)
stem(-M/2:M/2-1,fftshift(abs(Zk_phase(2,:))),'bs');
xlim([-20 20])
grid on
subplot(2,2,3)
stem(-M/2:M/2-1,fftshift(abs(Zk_phase(3,:))),'c*');
xlim([-20 20])
grid on
subplot(2,2,4)
stem(-M/2:M/2-1,fftshift(abs(Zk_phase(4,:))),'rv');
xlim([-20 20])
grid on

disp('4*(fs/N)/1024')
%%
% figure
% stem(kpha,fftshift(abs(Zk_RF_phase)),'m.-')
% xlim([-20 20])
% grid on

%%
% ADC 进行重采样
fs_adc = 512e6;
dsr = fs/fs_adc;
si_ds = downsample(si,dsr);

N_ds = N;
M_ds = length(si_ds)/N_ds; 

si_ds_arr = reshape(si_ds,[N_ds,M_ds]);
yds_arr = fft(si_ds_arr);
%%
kax_ds = -N_ds/2:N_ds/2-1;
figure(3)
stem(kax_ds,fftshift(abs(yds_arr(:,1))),'k.-')
xlim([0 500])
grid on
%%
ks_RF_ds = ks_RF*dsr;
yk_RF_seq_ds = yds_arr(ks_RF_ds+1,:); 
yk_RF_phase_ds = atan2(imag(yk_RF_seq_ds),real(yk_RF_seq_ds));
%%
slope_ds = 2*pi*dk*dsr;
theta_dk_ds = slope_ds*(0:M_ds-1);
figure
plot(unwrap(yk_RF_phase_ds(4,:)-theta_dk_ds),'k-')
hold on
plot(unwrap(yk_RF_phase_ds(4,:)),'r-')
hold on
plot((yk_RF_phase_ds(4,:)),'b-')
grid on
ylim([-5 5])

%% 2, 仅仅产生相位，模拟测试
clc;clear;
L = 4;
f = 4;
b = 3;

bi = 3;
bq = 5;
t= 0:L/1000:L-1;
theta = cos(2*pi*f*t);
theta_b = cos(2*pi*f*t)+b;

s = atan(tan(theta));
s1 = atan(sin(theta+bq)./cos(theta+bi));

tdata_DAC = atan2(sin(theta),cos(theta));
x1 = atan2(sin(theta_b),cos(theta_b));
%%
figure(1)
plot(theta,'k-');
hold on
plot(s1,'r-')
hold off
%%
figure(2)
plot(theta,'k-')
hold on
plot(x1,'r-')
hold off
%%
figure(3)
scatter(sin(theta_b),cos(theta_b),'m.')
hold on
scatter(sin(theta),cos(theta),'bh')
axis equal

%%
% DAC buffer数据测试
clear; clc;
fs = 4.096e9; %采样率
dsr = 8; % 下抽取率
c = 3e8; % 电磁波在传输线中的速度
% Ld = 30e-2; %传输线长度
% buffer = readmatrix('d:\TES_test_data\20211030_09_multil_withmodule_8190_buffer.txt', 'Delimiter','\n');

buffer = load('DAC_Buffer.mat');
tdata_DAC = buffer.buffer;

%%
N=2^10; % FFT长度为1024
M = length(tdata_DAC)/N; % buffer的FFT帧数

tdata_blk = reshape(tdata_DAC,[N,M]);
y_blk = fft(tdata_blk);
ks_RF = [2,32,36,40];
% p_dly = 2*pi*ks_RF*fs/N*Ld/c; % 不同频率的附加相位

tdata_ds_blk = reshape(downsample(tdata_DAC,8),[N,M/8]);
y_ds_blk = fft(tdata_ds_blk);

%%
kax = -N/2:N/2-1;
figure
stem(kax,fftshift(abs(y_blk(:,1))),'k.-')
xlim([0 50])
grid on

%%

y_phase = atan2(imag(y_blk(ks_RF+1,:)),real(y_blk(ks_RF+1,:)));
% y_ds_phase = atan2(imag(y_ds_blk(ks_RF*8+1,:)),real(y_ds_blk(ks_RF*8+1,:)));
figure
subplot(2,2,1)
plot(y_phase(1,:),'k.-')
% xlim([-20 20])
grid on
title('DAC-buffer 8MHz的调制相位波形')
subplot(2,2,2)
plot(y_phase(2,:),'b.-')
% xlim([-20 20])
grid on
title('DAC-buffer 128MHz的调制相位波形')
subplot(2,2,3)
plot(y_phase(3,:),'m.-')
% xlim([-20 20])
grid on
title('DAC-buffer 144MHz的调制相位波形')
subplot(2,2,4)
plot(y_phase(4,:),'r.-')
% xlim([-20 20])
grid on
title('DAC-buffer 160MHz的调制相位波形')
%%
Z_phase = fft(y_phase,[],2); % 沿行FFT，得到相位的变化
fM = fs/N/M*(-M/2:M/2-1)/1e3;
figure
subplot(2,2,1)
stem(fM,fftshift(abs(Z_phase(1,:))),'k.');
xlim([-50 50])
xlabel('kHz')
grid on
title('8MHz 的调制相位变化频率')
subplot(2,2,2)
stem(fM,fftshift(abs(Z_phase(2,:))),'bs');
xlim([-50 50])
xlabel('kHz')
grid on
title('128MHz 的调制相位变化频率')
subplot(2,2,3)
stem(fM,fftshift(abs(Z_phase(3,:))),'m*');
xlim([-50 50])
xlabel('kHz')
grid on
title('144MHz 的调制相位变化频率')
subplot(2,2,4)
stem(fM,fftshift(abs(Z_phase(4,:))),'rv');
xlim([-50 50])
xlabel('kHz')
grid on
title('160MHz 的调制相位变化频率')
%% 2021/11/9
% ADC 数据测试
clear;clc;
buffer = load('DAC_Buffer.mat');
tdata_DAC = buffer.buffer;

ADC = load('ADC_RCV_26214400.mat');
tdata_ADC = ADC.data;
clear ADC;

% tdata_DAC_downsampled = downsample(tdata_DAC,8);
% tdata_ADC_section = tdata_ADC(1:65536);

%%
% [c,lags] = xcorr(tdata_DAC_downsampled,tdata_ADC_section,'normalized');
% 
% figure
% stem(lags,c)
% xlim([-120 120])
%%
% figure
% plot(tdata_ADC_section,'b-')
% % hold on
% % plot(tdata_DAC_downsampled,'k.-')
% xlim([0 100])
%%
fs = 4.096e9; % RF信号的采样率
N = 1024; % FFT帧数
Mr = length(tdata_ADC)/N; % ADC 重采样后的FFT帧数
ks_RF = [2,32,36,40]; % DAC 发送的RF 频率 bin-index， 对应1024点FFT，4.096G/s采样
fs_RF = ks_RF/N*fs; % DAC 发送的RF 频率 (Hz)
tdata_rcv_blk = reshape(tdata_ADC, [N,Mr]);
y_rcv_blk = fft(tdata_rcv_blk);
dsr = 8;
ks_rcv_RF = ks_RF*dsr;
y_rcv_phase = atan2(imag(y_rcv_blk(ks_rcv_RF+1,:)),real(y_rcv_blk(ks_rcv_RF+1,:)));

%%
kax = -N/2:N/2-1;
figure
stem(kax,fftshift(abs(y_rcv_blk(:,1))),'k.-')
xlim([0 400])
grid on
%%
% unwrap调制相位，并去掉线性趋势
y_rcv_phase_detrend = zeros(4,Mr);
for ii=1:4
    y_rcv_phase_detrend(ii,:) = detrend(unwrap(y_rcv_phase(ii,:)));    
end
%%
% 画图比较去掉趋势项前后的调制相位
ksel = 1;
figure
subplot(2,1,1)
plot(y_rcv_phase_detrend(ksel,:),'r-');
hold on
plot(y_rcv_phase(ksel,:),'b-');
xlim([5500 6000])
title('在偏置点变化前后的调制相位波形')
legend('去掉线性趋势后的调制相位','原始数据调制相位')
grid on
subplot(2,1,2)
plot(y_rcv_phase_detrend(ksel,:),'r-');
hold on
plot(y_rcv_phase(ksel,:),'b-');
title('整体调制相位波形')
legend('去掉线性趋势后的调制相位','原始数据调制相位')
grid on
%
Nfft = N;
frcv = fs/dsr/N/Nfft*(-Nfft/2:Nfft/2-1)/1e3;
y_rcv_phase_seg = y_rcv_phase_detrend(ksel,1:Nfft);
Z_rcv_phase_seg = fft(y_rcv_phase_seg);
figure
plot(frcv,fftshift(abs(Z_rcv_phase_seg)),'k.-');
xlabel('kHz')
grid on

%%
%从数据文本文件中读取ADC数据
clear;clc;
filename = 'd:\TES_test_data\20211030_09_multil_withmodule_8190_shen_1s_cut_513828864.txt';
fd = fopen(filename);
formatSpec = '%d';

C_text = textscan(fd,formatSpec,Ns,'Delimiter','');
data = C_text{1,1};
fclose(fd);
save('ADC_RCV_26214400','data')