%%
% L-pnt x(n) sequence:
% signal frequency: integer+fraction  
% M-CoarseFFT & K-ZoomFFT phase variation and consistent correction.

clear;clc;
L = 20000; % original signal length.
M = 512; % CoarseFFT-length
K = 16; % ZoomFFT-length
n = 0:L-1;   % sampling index (time)
k = -M/2:M/2-1; % FFT-bin index (normalized frequency) 
srate = 1e9; % sampling rate.
kv = 32; % signal frequency integer-part for M-pnt FFT.
dkv = 0.24; % signal frequency fraction-part for M-pnt FFT.
% dkv = 0;
kvf = kv + dkv;

x = exp(1i*2*pi*kvf/M.*n); 
s = 1.5*exp(1i*2*pi*kv/M.*n); 

rect_win = (rectwin(M))'; % hamming window function

x_win = rect_win.*x(1:M);
X_win = fft(x_win);
X_win_mag = abs(X_win);

s_win = rect_win.*s(1:M);
S_win = fft(s_win);
S_win_mag = abs(S_win);

% plot time squence data.
figure(1),clf
subplot(2,1,1)
plot(0:M-1,real(x(1:M)),'k-','linew',2);
% ylim([0 3])
% xlim([500 1500])
grid on
title('Orig. signal sequence')

subplot(2,1,2)
plot(0:M-1,abs(x_win),'b-','linew',2);
grid on
title('Windowed signal sequence')
% ylim([0 3])
% xlim([500 1500])
%%
figure(2),clf
plot(k,fftshift(X_win_mag),'ks-','linew',1,'markerfacecolor','k','markersize',5);
hold on
stem(k,fftshift(S_win_mag),'bo-','linew',1,'markerfacecolor','b','markersize',5);
grid on
title('Signal spectrum')
xlim([0 100])
xlabel('FFT bin')
ylabel('FFT magnitude: (a.u.)')
legend({'\delta=0.24';'\delta=0'},'location','northwest')

%%
%
num_frames = floor(L/M);
x_arr = reshape(x(1:num_frames*M),[M,num_frames]);

Y_frames = fft(x_arr);

% Yp_mag_frames = abs(Yp_frames);

Y_kv_frames = Y_frames(kv+1,:);

m = 0:num_frames-1;

phase_rev = -2*pi*dkv.*m;
angles_kv = unwrap(angle(Y_kv_frames));

% Yp_kv_frames_rev = Yp_kv_frames.*(-1).^(-m*kv).*exp(-1i*pi*m*dkv);
angles_kv_rev = angles_kv + phase_rev;
%%
figure(3),clf
hold on
title(['Phase variations in M-FFT streamline: ', '; \delta=',num2str(dkv)])
plot(m,angles_kv,'ks-','linew',1,'markerfacecolor','k','markersize',5)
plot(m,angles_kv_rev,'bo-','linew',1,'markerfacecolor','b','markersize',5)
xlabel(['Sampling points at the',num2str(kv), '-th bin.(a.u.)'])
ylabel('Phase angle (rad)')
legend({'\delta=0.24';'\delta=0'},'location','northwest')
grid on
% ylim([0 2])
% plot(m,phasor_lags,'b-','linew',1)
%%


%% done
