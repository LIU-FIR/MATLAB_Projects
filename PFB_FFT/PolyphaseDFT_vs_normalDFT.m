%%
% 验证长点数的FFT和多相结构的FFT的等价
clear;clc;
N = 2048; % Signal sequence length
M = 256; % FFT-length
R = 8; % Taps in each polyphase branch. 
n = 0:N-1;   % sampling index (time)
k = -N/2:N/2-1; % FFT-bin index (normalized frequency) 
srate = 1e3; % sampling rate.
% f1 = 100; % signal frequency.
% f2 = 200;
qv1 = 120;
qv2 = 560;
% complex sinsoid signal (input)
x = exp(1i*2*pi*qv1/N.*n)+1.5*exp(1i*2*pi*qv2/N.*n); 
% win = (rectwin(N))'; 
ham_win = (hamming(N))'; % hamming window function
sinc_win = fir1(N-1,1/M,rectwin(N)); % sinc window function
sys_win = sinc_win.*ham_win; % composed system window function
% win = (gausswin(N))';

x_win = ham_win.*x;
%%
X_win = fft(x_win);
X_win_mag = abs(X_win);
% H = fft(ham_win);
% H_mag = abs(H);

figure(1),clf
subplot(2,1,1)
plot(n,real(x_win),'b-',n,imag(x_win),'r-');
hold on
plot(n,real(x),'b--',n,imag(x),'r--');
grid on
xlim([0 100])
title('Signal sequence')
subplot(2,1,2)
plot(k,fftshift(X_win_mag),'k-');
grid on
title('Signal spectrum')
% xlim([100 300])
% figure(1)
% subplot(2,1,1)
% plot(n,ham_win,'k-');
% grid on
% title('Window function')
% subplot(2,1,2)
% plot(k,fftshift(H_mag),'k-');
% grid on
% title('Window function spectrum')



%%
ham_win_pmat = reshape(ham_win,[M R]);% hamming window polyphase matrix
sys_win_pmat = reshape(sys_win,[M R]);% system window polyphase matrix
x_pmat = reshape(x,[M R]); % signal polyphase matrix
z = sum(ham_win_pmat.*(x_pmat),2);
Yp0 = fft(z); % use matrix-method to compute polyphase FFT.
Yp= pfft(x,ham_win,M,R); % use pfft() to compute polyphase FFT.
X_ds = downsample(X_win,R);% resample the orignal long FFT.

kd = (0:M-1);

%%
% plot original N-FFT spectrum and M-PFFT vs R-decimated N-PFFT
[~, idx_X_sorted] = sort(abs(X_win),'descend');
[~, idx_Xds_sorted] = sort(abs(X_ds),'descend');
figure(2),clf
% plot(kd,abs(Yp0),'ko')
% hold on
subplot(121)
plot(0:N-1,abs(X_win),'k-','markersize',10,'linew',1)
grid on
title([num2str(N),'-FFT magnitude spectrum'])
xlabel('FFT bins')
ylabel('FFT magnitude: (a.u.)')
set(gca,'xtick',[0 idx_X_sorted(2:-1:1)-1,N/2,N-1])
subplot(122)
hold on
plot(kd,abs(Yp),'rs','markersize',10,'linew',1)
plot(kd,abs(X_ds),'b*','markersize',8,'linew',1)
grid on
title([num2str(M),'-PFFT magnitude spectrum'])
set(gca,'xtick',[0 idx_Xds_sorted(2:-1:1)-1,M/2,M-1])
legend({[num2str(M),'-PFFT'];[num2str(R),'-decimated ',num2str(N),'-FFT']})
xlabel('PFFT bins')
ylabel('PFFT magnitude: (a.u.)')

%%
% PFB的频率响应：在某个bin内生成多个复正弦信号，以扫频方式画出每个正弦信号的频率响应。
kv = 16; % The selected FFT(M-pnts) bin
ks = kv-3:0.01:kv+3; % frequencies to be computed.
tN = (0:N-1); % original N-pnt length sample indices
tM = (0:M-1); % M-pnt length sample indices.

% Rp_frex: polyphase-FFT freq.response from the N-pnt sequence.  
% Rp_frex_win: windowed polyphase-FFT freq.response from the N-pnt sequence.
% R_frex: FFT freq.response from the M-pnt sequence.
% R_frex_win: windowed FFT freq.response from the M-pnt sequence.
[Rp_frex,Rp_frex_win,R_rex,R_frex_win] = deal(zeros(length(ks),1));
ham_win_Mpts = (hamming(M))';
for kk = 1:length(ks)
    xk = exp(1i*2*pi*ks(kk)/M*tN);   
    Yk = pfft(xk,ham_win,M,R); % use pfft()to compute M-channel from N-pts sequence.
    Yk1 = pfft(xk,sys_win,M,R);
    Rp_frex(kk)=abs(Yk(kv+1)); % pick up the kv-th bin. 
    Rp_frex_win(kk)=abs(Yk1(kv+1));
    x_sh = exp(1i*2*pi*ks(kk)/M*(0:M-1));
    x_win_sh = exp(1i*2*pi*ks(kk)/M*(0:M-1)).*ham_win_Mpts;
    Y_sh = fft(x_sh);% use fft()to compute M-channel from M-pts sequence.
    Y_win_sh = fft(x_win_sh);
    R_rex(kk) = abs(Y_sh(kv+1));% pick up the kv-th bin.
    R_frex_win(kk) = abs(Y_win_sh(kv+1));
end
disp('ends here')


%%
% 画出256点加窗FFT和2048点PFFT加窗
figure(3),clf
hold on
plot(ks, mag2db(R_rex/max(R_rex)),'k-','LineWidth',1)
plot(ks, mag2db(R_frex_win/max(R_frex_win)),'b-','LineWidth',1)
plot(ks, mag2db(Rp_frex/max(Rp_frex)),'m-','LineWidth',1)
plot(ks, mag2db(Rp_frex_win/max(Rp_frex_win)),'r-','LineWidth',1)
ylim([-80 1])
grid on
grid on
xlabel('FFT bins')
ylabel('dB')
% legend('256-pts FFT','256-pts FFT(hamming window)','2048-pts PFFT(hamming window)','2048-pts PFFT(hamming & sinc window)')

%%
sinc_win = fir1(N-1,1/M,rectwin(N));
figure(8)
plot(n,sinc_win,'b-')
% xk = exp(1i*2*pi*kv/M*t);
% Yk = polyphase_fft(xk,win,M,R);
% figure(5)
% plot(0:M-1,abs(Yk))
% xlim([0 50])