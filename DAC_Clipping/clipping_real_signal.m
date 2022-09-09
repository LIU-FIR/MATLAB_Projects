%%
clear;clc;
SteelBlue = '#36648B';
%% 
% 形成实载波信号
N = 1024;% 1帧FFT长度
n = 0:N-1;
fs = 2e9; % 采样率
Nc = N/16; % 载波数量

% ks = [20 64 100 120]; % 载波对应的bins
ks = randperm(N/2-1,Nc);% 载波对应的bins
phases = randperm(N,Nc)/N*2*pi; % 从 [1,N]、N*2*pi中选出Nc个不同的相位
% phases = randperm(N,Nc)/N*2*pi-pi; % 从 [1,N]、N*2*pi中选出Nc个不同的相位


AF = 2; % DAC full-scale
P_avg = AF.^2/2; % 单路ADC采样信号总功率
A = AF/Nc^0.5; % r.m.s值平方为A^2/2，故计算幅度和复信号相同

%%
ind_ks = ks+1;

% 产生不同相位和频率的载波信号
s_arr = zeros(N,Nc);

for vv = 1:Nc
    s_arr(:,vv) = A*cos(2*pi*ks(vv)*n/N + phases(vv)); 
end
Xs = 1/N*fft(s_arr);

s_sum = sum(s_arr,2);
X_sum = 1/N*fft(s_sum);

s_DAC_clip = s_sum;
s_DAC_clip(s_DAC_clip>AF) = AF;
s_DAC_clip(s_DAC_clip<-AF) = -AF;

PX_sum = abs(X_sum).^2;

X_DAC_clip = 1/N*fft(s_DAC_clip);
PX_DAC_clip = abs(X_DAC_clip).^2; 

PX_sum(PX_sum==0)=eps;
PX_DAC_clip(PX_DAC_clip==0)=eps;


%%
figure(1)
subplot(2,1,1)
plot(s_sum,'k-','linew',1)
hold on
plot([1,N],AF*[1,1],'r--','linew',1)
plot([1,N],AF*[-1,-1],'r--','linew',1)
title('The sum of all carriers (in time)')
% xlim([0 N/2])
caxis = gca;
subplot(2,1,2)
hold on
plot(n,s_DAC_clip,'b-','linew',1)
plot([1,N],AF*[1,1],'r--','linew',1)
plot([1,N],AF*[-1,-1],'r--','linew',1)
title('DAC-clipping of all carriers'' sum')
% xlim([0 N/2])
ylim(caxis.YLim)
%
figure(2)
subplot(2,1,1)
plot(pow2db(PX_sum/P_avg),'k-','linew',1)
hold on
plot(ind_ks,pow2db(PX_sum(ind_ks)/P_avg),'r.')
xlim([0,N/2])
xlabel('FFT bins [0,N/2]')
ylabel('dB')
title('The power spectrum of all carriers'' sum')
subplot(2,1,2)
plot(pow2db(PX_DAC_clip/P_avg),'b-','linew',1)
hold on
plot(ind_ks,pow2db(PX_DAC_clip(ind_ks)/P_avg),'r.')
xlim([0,N/2])
xlabel('FFT bins [0,N/2]')
ylabel('dB')
title('The DAC-clipping power spectrum of all carriers'' sum')
%%
% 构造实陷波滤波器
p = zeros(N,1);
p(1) = 1;
X_p = fft(p);
MX_p = abs(X_p);

% p_tmp = ifft(X_pr);
X_pr = X_p;
X_pr(ind_ks) = 0;
X_pr(N-ind_ks+1) = 0; % 构造实滤波脉冲时需要剔除关于bin=N/2对称的频率分量

MX_pr = abs(X_pr);
pr = ifft(X_pr);
%%
figure(3)
subplot(2,2,1)
plot(p,'k-','linew',1)
subplot(2,2,3)
plot(MX_p,'k-','linew',1)
ylim([0 2])
subplot(2,2,4)
plot(MX_pr,'b-','linew',1)
subplot(2,2,2)
plot(real(pr),'b-','linew',1)
hold on
plot(imag(pr),'r-','linew',1)

%%
s_AG_clip = s_sum;
ft_pr = fft(real(pr),2*N);
n_iter = 0;

tol = 1e6*eps;

while(max(s_AG_clip)-AF >tol  || min(s_AG_clip)+AF< -tol)
    c = zeros(N,1);
    ind_upper = s_AG_clip > AF;
    ind_lower = s_AG_clip < -AF;
    c(ind_upper) = AF - s_AG_clip(ind_upper);
    c(ind_lower) = -AF - s_AG_clip(ind_lower);

    ft_c = fft(c,2*N);
    ft_cc = ft_c.* ft_pr;   

    cc_full = real(ifft(ft_cc));
    cc = cc_full(1:N);
    % cc = conv(c,pulse_t_rev,'same');
    s_AG_clip = s_AG_clip + cc;
    n_iter = n_iter+1;
    disp(n_iter);
end
%%
X_AG_clip = 1/N*fft(s_AG_clip);
PX_AG_clip = abs(X_AG_clip).^2;

phases_AG_clip = angle(X_AG_clip(ind_ks))'; % 计算AG算法 clipping后的载波频点相位
phases_AG_clip = mod(phases_AG_clip,2*pi); % 相位取到[0, 2*pi]内

phases_DAC_clip = angle(X_DAC_clip(ind_ks))'; % 计算DAC clipping后的载波频点相位
phases_DAC_clip = mod(phases_DAC_clip,2*pi); % 相位取到[0, 2*pi]内

err1_phase_ratio = (phases_DAC_clip - phases)./2/pi;
err1_phase_abs = phases_DAC_clip - phases;
% plot(err1_phase)

err2_phase_ratio = (phases_AG_clip - phases)./2/pi;
err2_phase_abs = phases_AG_clip - phases;

% disp(err1_phase)
% disp(err2_phase)
%%
figure(4)
subplot(2,1,1)
hold on
plot(s_sum,'k-','linew',1)
plot([1,N],AF*[1,1],'r--','linew',1)
plot([1,N],AF*[-1,-1],'r--','linew',1)
title('The sum of all carriers (in time)')
caxis = gca;
subplot(2,1,2)
hold on
phd3 = plot(s_AG_clip,'-','linew',1);
phd3.Color = SteelBlue;
plot([1,N],AF*[1,1],'r--','linew',1)
plot([1,N],AF*[-1,-1],'r--','linew',1)
title('AG-clipping of all carriers'' sum')
% xlim([0 N/2])
ylim(caxis.YLim)
%
figure(5)
subplot(2,1,1)
plot(pow2db(PX_sum/P_avg),'k-','linew',1)
hold on
plot(ind_ks,pow2db(PX_sum(ind_ks)/P_avg),'r.')
xlim([0 N/2])
xlabel('FFT bins [0,N/2]')
ylabel('dB')
title('The power spectrum of all carriers'' sum')
% subplot(2,2,3)
% plot(angle(X_sum),'k-','linew',1)
% xlim([0 N/2])
subplot(2,1,2)
phdl = plot(pow2db(PX_AG_clip/P_avg),'-','linew',1);
% phdl.Color = '#696969';
phdl.Color = SteelBlue;
hold on
plot(ind_ks,pow2db(PX_AG_clip(ind_ks)/P_avg),'r.')
xlim([0 N/2])
xlabel('FFT bins [0,N/2]')
ylabel('dB')
title('The AG-clipping power spectrum of all carriers'' sum')
% subplot(2,2,4)
% plot(angle(X_AG_clip),'b-','linew',1)
% xlim([0 N/2])
%%
figure(6)
subplot(2,1,1)
plot(err1_phase_ratio,'b.-','linew',1);
caxis = gca;
% plot(err1_phase_abs,'b.-','linew',1)
title('The relative phase error between DAC-clipping carriers and unclipped carriers')
% title('The absolute phase error between DAC-clipping carriers and unclipped carriers')
% hold on
subplot(2,1,2)
phd2 = plot(err2_phase_ratio,'.-','linew',1);
% phd2 = plot(err2_phase_abs,'.-','linew',1);
phd2.Color = SteelBlue;
ylim(caxis.YLim)
title('The relative phase error between AG-clipping carriers and unclipped carriers')
% title('The absolute phase error between AG-clipping carriers and unclipped carriers')

%%
figure(7)

phd1 = plot(err1_phase_ratio,'b.-','linew',1);
phd1.MarkerSize=15;
hold on
% plot(err1_phase_abs,'b.-','linew',1)
phd2 = plot(err2_phase_ratio,'.--','linew',1);
% phd2 = plot(err2_phase_abs,'.-','linew',1);
phd2.Color = SteelBlue;
phd2.MarkerSize=15;
title('The relative phase error between clipped and unclipped carriers')
legend('DAC-clipping','AG-clipping')
% title('The absolute phase error between AG-clipping carriers and unclipped carriers')

