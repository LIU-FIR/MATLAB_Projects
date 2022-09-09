%%
% 验证分数频率的DFT在相近bin上输出特性
% 将长时间序列划分为多个N点序列帧，

clear;clc;
N = 512;  % 第一级FFT长度
W = 100;  % FFT帧数
n = 0:N*W-1;   % 信号持续长度
k = -N/2:N/2-1;
kv = 200.25; % N点信号上的频率索引，非整数
kvi = 200;
dkv = kv-kvi;
phi0 = pi/8; % 引入初始相位
phi1 = pi/16;
phi_add = pi*dkv*(N-1)/N;
x = exp(1i*2*pi*kvi/N*n+1i*phi0) + exp(1i*2*pi*kv/N*n+1i*phi1);
% x = exp(1i*2*pi*kv/N*n+1i*phi1);
xp = reshape(x,[N,W]);
%%
yp = fft(xp);
ym_kv = yp(kvi+1,:);
m = 0:W-1; % 不同帧形成的时序

figure(1)
subplot(2,1,1)
plot(k,fftshift(abs(yp(:,1))),'k-o');
grid on
xlim([180 220])
title("k1=200和k2=200.25两个频率的信号频谱")
% title('Window function')
subplot(2,1,2)
plot(m,real(ym_kv),'b-',m,imag(ym_kv),'r-');
title("100帧FFT的连续输出序列，在第200 bin上观察")

legend('real','imag')
grid on


%%
L = 128;
z = fft(ym_kv,L);
r = -L/2:1:L/2-1;

%%

phase_rotator = exp(-1i*2*pi*dkv*m); %抵消差频在100帧上产生的相位变化
s = ym_kv.*phase_rotator;
S = fft(s,L);
figure(3)
subplot(2,1,1)
plot(r,fftshift(abs(z)),'k')
grid on
title("将给定bin的100帧序列直接做FFT的结果")
% figure(4)
subplot(2,1,2)
plot(r,fftshift(abs(S)),'b')
title("将给定bin的100帧序列相位停止后做FFT的结果")
grid on
%%
% phi_c = atan2(imag(z(32)),real(z(32)));
phi_c = atan2(imag(S(1)/L),real(S(1)/L));
%%
q = 0:2*pi/8:2*pi;
a = exp(1i*q);
figure(4)

for i = 1:10
    scatter(imag(ym_kv)./abs(ym_kv),real(ym_kv)./abs(ym_kv),'b')
    hold on
    scatter(imag(ym_kv(i))./abs(ym_kv(i)),real(ym_kv(i))./abs(ym_kv(i)),'filled','r')
    rectangle('Position',[-1, -1, 2, 2],'Curvature',[1, 1]);axis equal;
    hold off
    pause(0.5)
end
