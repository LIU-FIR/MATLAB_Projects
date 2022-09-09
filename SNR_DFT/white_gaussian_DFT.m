
%%
% 计算高斯白噪声在FFT之后每个bin上的功率
% 利用多个独立同分布的高斯白噪声数据的FFT后功率谱进行估计
clc;clear;

M = 10000;
Nfft = 4096;
L = M*Nfft;
sigma = 5;
wgnoise = sigma.*randn(Nfft,M); 

%
X = fft(wgnoise)/Nfft;
P = real(X).^2 + imag(X).^2;

P_estimate = mean(P,2); % 利用不同数据窗的频谱均值估计P
%

mean_P_sta = sigma^2/Nfft; % 计算P的均值统计值，服从Chi-square分布
var_P_sta = 2*sigma^2/Nfft; % 计算P的方差统计值，服从Chi-square分布

% error_ratio = abs((mean_P_sta - P_estimate)/P_estimate);
% disp([mean_P_sta, P_estimate, error_ratio])
%
fk = -Nfft/2:Nfft/2-1;
figure
plot(fk,fftshift(P_estimate),'k.','linew',2)
hold on
% plot([fk(1) fk(end)],[1 1]*P_estimate,'r-','linew',1.5)
% hold on
plot([fk(1) fk(end)],[1 1]*mean_P_sta,'b-','linew',1.5)
ylim([mean_P_sta-var_P_sta/4 mean_P_sta+var_P_sta/4])

grid on
xlabel('DFT Bins: ','FontSize',14','Interpreter','latex')
ylabel('$\mathrm{Power: |W(k)|^2}$','Interpreter','latex','FontSize',14)
title({'Power Spectrum Density Estimation'},{ '$\mathrm{Gaussian\, noise\, sequence\,\sigma = 5, L = 4096{\times}10000, N_{FFT} = 4096}$'},'Interpreter','latex','FontSize',14)
legend('Mean Power estimation over all bins','Mean Power (from Chi-square distribution)','FontSize',10,'Interpreter','latex')
% ylim([mean_P_sta-error_ratio*mean_P_sta*5 mean_P_sta+error_ratio*mean_P_sta])
%%
% 计算高斯白噪声在FFT之后每个bin上的功率
% 利用单个独立同分布的高斯白噪声数据窗，假设序列为平稳过程
clc;clear;

M = 1;
Nfft = 4096;
L = M*Nfft;
sigma = 5;
wgnoise = sigma.*randn(Nfft,M); 

%
X = fft(wgnoise)/Nfft;
P = real(X).^2 + imag(X).^2;

mean_P_sta = sigma^2/Nfft; % 计算P的均值统计值，服从Chi-square分布
var_P_sta = 2*sigma^2/Nfft; % 计算P的方差统计值，服从Chi-square分布
%
X = fft(wgnoise)/Nfft;
P = real(X).^2 + imag(X).^2;

P_estimate = mean(P); % 利用所有bin的功率谱估计P
%

mean_P_sta = sigma^2/Nfft; % 计算P的均值统计值，服从Chi-square分布
var_P_sta = 2*sigma^2/Nfft; % 计算P的方差统计值，服从Chi-square分布

error_ratio = abs((mean_P_sta - P_estimate)/P_estimate);
disp([mean_P_sta, P_estimate, error_ratio])
%
fk = -Nfft/2:Nfft/2-1;
figure
subplot(2,1,1)
plot(fk,fftshift(P),'k.','linew',2)
hold on
plot([fk(1) fk(end)],[1 1]*P_estimate,'r-','linew',1.5)
hold on
plot([fk(1) fk(end)],[1 1]*mean_P_sta,'b-','linew',1.5)
% ylim([mean_P_sta-var_P_sta/2 mean_P_sta+var_P_sta/2])

grid on
xlabel('DFT Bins: ','FontSize',14','Interpreter','latex')
ylabel('$\mathrm{Power: |W(k)|^2}$','Interpreter','latex','FontSize',14)
title({'Power Spectrum Density Estimation'},{ '$\mathrm{Gaussian\, noise\, sequence\,\sigma = 5, L = 4096, N_{FFT} = 4096}$'},'Interpreter','latex','FontSize',14)
legend('$|W(k)|^2$ on all bins','Mean Power estimation over all bins','Mean Power (from Chi-square distribution)','FontSize',10,'Interpreter','latex')

subplot(2,1,2)
plot(fk,fftshift(P),'k.','linew',2)
hold on
plot([fk(1) fk(end)],[1 1]*P_estimate,'r-','linew',1.5)
hold on
plot([fk(1) fk(end)],[1 1]*mean_P_sta,'b-','linew',1.5)
ylim([mean_P_sta-error_ratio*mean_P_sta*5 mean_P_sta+error_ratio*mean_P_sta*5])
xlabel('DFT Bins: ','FontSize',14','Interpreter','latex')
ylabel('$\mathrm{Power: |W(k)|^2}$','Interpreter','latex','FontSize',14)
legend('$|W(k)|^2$ on all bins','Mean Power estimation over all bins','Mean Power (from Chi-square distribution)','FontSize',10,'Interpreter','latex')
grid on
