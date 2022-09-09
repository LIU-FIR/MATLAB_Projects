%%
% 验证被信息信号(message signal)调制的信号特性
% 

clear;clc;
fs = 200;
fm = 10;
fc = 50;
t = 0:1/fs:1;
x = sin(2*pi*fm*t) + randn(size(t))/10;

y = modulate(x,fc,fs,'am');
N = 256;
X = fft(x.*hamming(length(t))',N);
Y = fft(y.*hamming(length(t))',N);
k = -N/2:N/2-1;
Xmag = abs(X);
Ymag = abs(Y);
% pwelch([x;y]',hamming(100),80,1024,fs,'centered')

%%
figure(1)
subplot(2,1,1)
plot(t,x,'k-')
grid on
subplot(2,1,2)
plot(t,y,'b-')
grid on
%%
figure(2)
subplot(2,1,1)
% plot(k,fftshift(Xmag),'k-')
plot(k/N*fs,fftshift(Xmag),'k-')
grid on
subplot(2,1,2)
% plot(k,fftshift(Ymag),'b-')
plot(k/N*fs,fftshift(Ymag),'b-')
grid on