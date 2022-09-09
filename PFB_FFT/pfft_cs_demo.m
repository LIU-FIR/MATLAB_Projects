%%
% 验证长序列FFT多相结构的连续输出
clear;clc;

M = 8192; % Long-FFT-length
R = 4; % Taps in each polyphase branch. 
N = M/R; % Polyphase FFT-length
L = M*120; % original signal length.
n = 0:L-1;   % sampling index (time)
P = 64; % 2-nd FFT-length
OSR = 1; % over-sampled ratio

ksM = -M/2:M/2-1; % FFT-bin index (normalized frequency) 
ksP = -P/2:P/2-1; % FFT-bin index (normalized frequency) 
srate = 2048e6; % sampling rate.

kM = 160; % signal frequency central-bin for M-pnt long-FFT. kM in [-M/2:M/2-1]
kP = 31; % signal frequency central-bin in second-FFT. kP in [-P/2:P/2-1]
% dkM = kP*R*OSR/P ;% signal frequency fraction-part for M-pnt long-FFT. kP in [-M/2:M/2-1]
dkM = kP*OSR/P ;% signal frequency fraction-part for M-pnt long-FFT. kP in [-M/2:M/2-1]

kvf = kM + dkM; % signal real frequency for M-pnt FFT.
fv = srate/M*(kM + kP/P*OSR);
dfv = srate/M*(1 + 1/P*OSR);

disp(['Sampling Rate: ',num2str(srate),' /s'])
disp(['kM-bin (M-FFT): ',num2str(kM),', kP-bin (2nd-P-FFT): ',num2str(kP)])
disp(['Real-frequency-bin(M-FFT): ',num2str(kvf)])
disp(['Real-frequency: ',num2str(fv), ' Hz'])
disp(['frequency resolution: ',num2str(dfv), ' Hz'])

%%
% x = zeros(1,L);
% for vv = 1:P
%     x = exp(1i*2*pi*kvf(vv)/M.*n) + x; 
% end
x = exp(1i*2*pi*fv/srate.*n); 
% win = (rectwin(N))'; 
rect_win = ones(1,M)./M;
ham_win = (hamming(M))'; % hamming window function
ham_win = ham_win./sum(ham_win);
% sinc_win = fir1(N-1,1/M,rectwin(N)); % sinc window function

% sinc_win = sinc(pi.*((0:M-1)-M/2)./N);
% sinc_win = sinc(pi.*((0:M-1)-M/2)./(2*N))/2;
sinc_win = sinc(pi.*((0:M-1)-M/2)./(4*N))/4;
sinc_win = sinc_win./sum(sinc_win);
%
sys_win = sinc_win.*ham_win; % composed system window function
sys_win = sys_win./sum(sys_win);
sys_win1 = rect_win;

%
% plot M-pnt windowed time squence data .
x_sysWin = x(1:M).* sys_win;
%%
figure(1),clf
% subplot(211)
hold on
plot(0:M-1,real(x_sysWin),'b-','linew',2);
% grid on
title('Windowed-carrier-signal sequence''s realpart')
% subplot(212)
plot(0:M-1,imag(x_sysWin),'r-','linew',2);
grid on

xlim([3000 5000])

%%


os_ratio = OSR; % the oversampled ratio. os_ratio>1 means that outputs from pfft_cs's bins refresh faster.
slide_len = M/os_ratio; % when being loaded,the start position in orig.sequence increase by slide_len

load_times = floor((L-M)/slide_len)+1; % deducted from geometric relation
disp(['load_times: ' num2str(load_times)])
%
Y_PFB_frames = zeros(N,load_times);
Y1_PFB_frames = zeros(N,load_times);

for cnt = 1:load_times
    % for every sliding, compute the current head and tail of the
    % ready-to-pfft_cs sequence.
    head_ptr = 1+(cnt-1)*slide_len;
    tail_ptr = 1+(cnt-1)*slide_len + M-1;

    xk_reg = x(head_ptr:tail_ptr); % load the ready-to-pfft_cs sequence in register.
    Y_PFB_frames(:,cnt) = pfft_cs(xk_reg,sys_win,N,R);
    Y1_PFB_frames(:,cnt) = pfft_cs(xk_reg,rect_win,N,R);
end    

disp('ends here')
%
% run 2-nd FFT to extract the real frequency
chan_Nidx = kM/R+1; %% Caution: N channel decimation from M-FFT
Y_PFB_samples = Y_PFB_frames(chan_Nidx,:);
Y1_PFB_samples = Y1_PFB_frames(chan_Nidx,:);

Z = fft(Y_PFB_samples(1:P));
Z1 = fft(Y1_PFB_samples(1:P));
%

%%
ksN = -N/2:N/2-1;
figure(3),clf
hold on
plot(ksN, fftshift(mag2db(abs(Y_PFB_frames(:,1)))),'k-','LineWidth',2)

figure(4),clf
hold on
plot(ksP,fftshift(mag2db(abs(Z)./max(abs(Z)))),'bs','LineWidth',2)
plot(ksP,fftshift(mag2db(abs(Z1)./max(abs(Z1)))),'ks','LineWidth',2)
ylim([-9 2])
% plot(ksP,(mag2db(abs(Z))),'bs','LineWidth',2)
% plot(mag2db(abs(Z)),'bs','LineWidth',2)
grid on
% ylim([-100 5])
%%
Y_PFB = pfft(x,sys_win,N,R,2);
Y_samples = Y_PFB(chan_Nidx,:);
Z0 = fft(Y_samples(1:P));
%%
figure(5),clf
hold on
plot(ksP,fftshift(mag2db(abs(Z)./max(abs(Z)))),'b.','LineWidth',2)
plot(ksP,fftshift(mag2db(abs(Z0)./max(abs(Z1)))),'ko','LineWidth',2)