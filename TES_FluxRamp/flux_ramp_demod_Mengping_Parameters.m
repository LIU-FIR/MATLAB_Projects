%%
% 根据梦萍的参数设置，给定DAC存储信号(fs=4096M/s)，ADC按照fs/8采样
% 
%% 2022-05-16
clc;clear;
Nx=2^10;Ny=2^9;
N = Nx*Ny; % DAC Buffer-length
n = 0:N-1;
srate_DAC = 4.096e9; % DAC sampling rate
srate_ADC = 512e6; % ADC sampling rate
Rd = srate_DAC/srate_ADC; % downsampled ratio

f1 = 8e6;
f2 = 128e6;
f3 = 144e6;
f4 = 160e6;
fc = 15.625e3; % modulate phase's frequency
A = 8190/4;

phi_mod = abs(sin(2*pi*fc/srate_DAC*n)); % modulated phase by fc-sinsoid carrier.
x = A*cos(2*pi*f1/srate_DAC*n + phi_mod)...
    + A*cos(2*pi*f2/srate_DAC*n)...
    + A*cos(2*pi*f3/srate_DAC*n - phi_mod)...
    + A*cos(2*pi*f4/srate_DAC*n);

%%
y = downsample(x,Rd);
m = 0:N/Rd-1;

figure(1),clf
subplot(211)
plot(n,x,'b-')
title('Signal in DAC buffer')
xlim([0 3200])
subplot(212)
plot(m,y,'k-')
title('Signal downsampled by ADC')
xlim([0 400])
%%
Nfft = 512;
y_win = y(1:Nfft);
Y = abs(fft(y_win));
k = (-Nfft/2:Nfft/2-1);
fk = (-Nfft/2:Nfft/2-1)*srate_ADC/Nfft/1e6; % MHz

figure(2),clf
plot(k,fftshift(Y),'k-','linew',1)
grid on
title('DFT magnitude of down-sampled signal')
xlim([0 Nfft/2-1])
%%
Ly = length(y); 
os_ratio = 1;
slide_len = Nfft/os_ratio; % sliding length to form a new frame.

N_frames = floor((Ly-Nfft)/slide_len)+1; % compute FFT-frame number of the down-sampled sequence  
Y_frames = zeros(Nfft,N_frames);

idx_frex = [8+1 128+1 144+1 160+1];

for cnt = 1:N_frames
    % for every sliding, compute the current head and tail of the
    % ready-to-pfft sequence.
    head_ptr = 1+(cnt-1)*slide_len;
    tail_ptr = 1+(cnt-1)*slide_len + Nfft-1;
    
    y_reg = y(head_ptr:tail_ptr); % load the ready-to-pfft sequence in register.
    Y_frames(:,cnt) = fft(y_reg);
end

Y_f1 = Y_frames(idx_frex(1),:); % 8MHz-bin, phase-modulated
Y_f2 = Y_frames(idx_frex(2),:);
Y_f3 = Y_frames(idx_frex(3),:); % 144MHz-bin, phase-modulated
Y_f4 = Y_frames(idx_frex(4),:);

disp(['N_frames: ',num2str(N_frames)])
%
figure(3),clf
subplot(221),hold on
plot(angle(Y_f1),'ks-','linew',1)
grid on
title('Phase at 8MHz-bin')
subplot(222),hold on
plot(angle(Y_f2),'bs-','linew',1)
ylim([-1 1])
title('Phase at 120MHz-bin')
grid on
subplot(223),hold on
plot(angle(Y_f3),'rs-','linew',1)
title('Phase at 144MHz-bin')
grid on
subplot(224),hold on
plot(angle(Y_f4),'ms-','linew',1)
title('Phase at 160MHz-bin')
ylim([-1 1])
grid on

