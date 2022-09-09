%%
% 验证长点数N-FFT等价为多个并行M-PFFT (N=M*R)
%
clear;clc;
N = 6144; % Signal sequence length
M = 1024; % FFT-length
R = 6; % Taps in each polyphase branch. 
n = 0:N-1;   % sampling index (time)

srate = 1e3; % sampling rate.

% signal frequency bins, which are even-indexed and odd-indexed
v1 = 300; % R*p 
v2 = 601; % R*p+1
v3 = 902; % R*p+2
v4 = 1203; % R*p+3
v5 = 1804; % R*p+4
v6 = 2405; % R*p+5

% complex sinsoid signal (input)

x = 0.6*exp(1i*2*pi*v1/N.*n)+1*exp(1i*2*pi*v2/N.*n)...
    + 1.5*exp(1i*2*pi*v3/N.*n)+2*exp(1i*2*pi*v4/N.*n)...
    + 2.5*exp(1i*2*pi*v5/N.*n)+3*exp(1i*2*pi*v6/N.*n);

ham_win = (hamming(N))'; % hamming window function
sinc_win = fir1(N-1,1/M,rectwin(N)); % sinc window function
sys_win = sinc_win.*ham_win; % composed system window function
sys_win = (rectwin(N))';

x_win = sys_win.*x;
X_win = fft(x_win);
X_win_mag = abs(X_win);

%%
k = -N/2:N/2-1; % FFT-bin index (normalized frequency) 
figure(1),clf
plot(k,fftshift(X_win_mag),'k-','linew',1);
grid on
title('Signal spectrum')
xlabel('FFT bins')
ylabel('Magnitude: (a.u.)')
xlim([0 N/2-1])

%%

X_ds0 = downsample(X_win,R);% resample the orignal long FFT.
X_ds1 = downsample(X_win,R,1);
X_ds2 = downsample(X_win,R,2);
X_ds3 = downsample(X_win,R,3);
X_ds4 = downsample(X_win,R,4);
X_ds5 = downsample(X_win,R,5);

% use pfft_n() to compute polyphase FFT.
Yp0 = pfft_n(x,sys_win,M,R,0); 
Yp1 = pfft_n(x,sys_win,M,R,1);
Yp2 = pfft_n(x,sys_win,M,R,2);
Yp3 = pfft_n(x,sys_win,M,R,3);
Yp4 = pfft_n(x,sys_win,M,R,4);
Yp5 = pfft_n(x,sys_win,M,R,5);

kd = (0:M-1);
%%

% figure(2),clf
% subplot(221),hold on
% plot(kd,abs(X_ds0),'b*-','markersize',8,'linew',1)
% grid on
% title(['Decimated FFT:',' DSI=',num2str(0)])
% xlabel('Decimated FFT bins')
% ylabel('FFT magnitude: (a.u.)')
% ylim([0 5000])
% 
% subplot(222),hold on
% plot(kd,abs(X_ds1),'b*-','markersize',8,'linew',1)
% grid on
% title(['Decimated FFT:',' DSI=',num2str(1)])
% xlabel('Decimated FFT bins')
% ylabel('FFT magnitude: (a.u.)')
% ylim([0 5000])
% 
% subplot(223),hold on
% plot(kd,abs(X_ds2),'b*-','markersize',8,'linew',1)
% grid on
% title(['Decimated FFT:',' DSI=',num2str(2)])
% xlabel('Decimated FFT bins')
% ylabel('FFT magnitude: (a.u.)')
% ylim([0 5000])
% 
% subplot(224),hold on
% plot(kd,abs(X_ds3),'b*-','markersize',8,'linew',1)
% grid on
% title(['Decimated FFT:',' DSI=',num2str(3)])
% xlabel('Decimated FFT bins')
% ylabel('FFT magnitude: (a.u.)')
% ylim([0 5000])
%%
figure(3),clf
subplot(231),hold on
plot(kd,abs(Yp0),'rs','markersize',10,'linew',1)
plot(kd,abs(X_ds0),'b*-','markersize',8,'linew',1)
grid on
title(['PFFT and decimated FFT:',' DSI=',num2str(0)])
% set(gca,'xtick',[0 idx_Xds_sorted(2:-1:1)-1,M/2,M-1])
legend({[num2str(M),'-PFFT'];[num2str(R),'-decimated ',num2str(N),'-FFT']})
xlabel('PFFT bins')
ylabel('PFFT magnitude: (a.u.)')
ylim([0 20000])

subplot(232),hold on
plot(kd,abs(Yp1),'rs','markersize',10,'linew',1)
plot(kd,abs(X_ds1),'b*-','markersize',8,'linew',1)
grid on
title(['PFFT and decimated FFT:',' DSI=',num2str(1)])
% set(gca,'xtick',[0 idx_Xds_sorted(2:-1:1)-1,M/2,M-1])
legend({[num2str(M),'-PFFT'];[num2str(R),'-decimated ',num2str(N),'-FFT']})
xlabel('PFFT bins')
ylabel('PFFT magnitude: (a.u.)')
ylim([0 20000])

subplot(233),hold on
plot(kd,abs(Yp2),'rs','markersize',10,'linew',1)
plot(kd,abs(X_ds2),'b*-','markersize',8,'linew',1)
grid on
title(['PFFT and decimated FFT:',' DSI=',num2str(2)])
% set(gca,'xtick',[0 idx_Xds_sorted(2:-1:1)-1,M/2,M-1])
legend({[num2str(M),'-PFFT'];[num2str(R),'-decimated ',num2str(N),'-FFT']})
xlabel('PFFT bins')
ylabel('PFFT magnitude: (a.u.)')
ylim([0 20000])

subplot(234),hold on
plot(kd,abs(Yp3),'rs','markersize',10,'linew',1)
plot(kd,abs(X_ds3),'b*-','markersize',8,'linew',1)
grid on
title(['PFFT and decimated FFT:',' DSI=',num2str(3)])
% set(gca,'xtick',[0 idx_Xds_sorted(2:-1:1)-1,M/2,M-1])
legend({[num2str(M),'-PFFT'];[num2str(R),'-decimated ',num2str(N),'-FFT']})
xlabel('PFFT bins')
ylabel('PFFT magnitude: (a.u.)')
ylim([0 20000])

subplot(235),hold on
plot(kd,abs(Yp4),'rs','markersize',10,'linew',1)
plot(kd,abs(X_ds4),'b*-','markersize',8,'linew',1)
grid on
title(['PFFT and decimated FFT:',' DSI=',num2str(4)])
% set(gca,'xtick',[0 idx_Xds_sorted(2:-1:1)-1,M/2,M-1])
legend({[num2str(M),'-PFFT'];[num2str(R),'-decimated ',num2str(N),'-FFT']})
xlabel('PFFT bins')
ylabel('PFFT magnitude: (a.u.)')
ylim([0 20000])

subplot(236),hold on
plot(kd,abs(Yp5),'rs','markersize',10,'linew',1)
plot(kd,abs(X_ds5),'b*-','markersize',8,'linew',1)
grid on
title(['PFFT and decimated FFT:',' DSI=',num2str(5)])
% set(gca,'xtick',[0 idx_Xds_sorted(2:-1:1)-1,M/2,M-1])
legend({[num2str(M),'-PFFT'];[num2str(R),'-decimated ',num2str(N),'-FFT']})
xlabel('PFFT bins')
ylabel('PFFT magnitude: (a.u.)')
ylim([0 20000])
%%
