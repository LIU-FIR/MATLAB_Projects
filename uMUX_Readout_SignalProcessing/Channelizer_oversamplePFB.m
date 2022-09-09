%%
% 验证长序列FFT多相结构的连续输出
clear;clc;

M = 8192; % Long-FFT-length
R = 4; % Taps in each polyphase branch. 
N = M/R; % Polyphase FFT-length
L = M*120; % original signal length.
n = 0:L-1;   % sampling index (time)
P = 64; % 2-nd FFT-length



srate = 2048e6; % sampling rate.

kM = 160; % signal frequency central-bin for M-pnt long-FFT. kM in [-M/2:M/2-1]
kP_critical = P/2-1; % signal frequency central-bin in second-FFT. kP in [-P/2:P/2-1]

% dkM = kP_critical/P ;% signal frequency fraction-part for M-pnt long-FFT. kP in [-M/2:M/2-1]

% kvf = kM + dkM; % signal real frequency for M-pnt FFT.
fkM_bin = srate/M*kM;
fv_critical = srate/M*(kM + kP_critical/P);
dfv_critical = srate/M*(1 + 1/P);

kP_sig = kP_critical+1;
fv_sig = srate/M*(kM + kP_sig/P);
% fv_sig = srate/M*kM;

disp(['Sampling Rate: ',num2str(srate),' /s'])
disp(['Bin-center frequency: ',num2str(fkM_bin),' /s'])
% disp(['kM-bin (M-FFT): ',num2str(kM),', kP-bin (2nd-P-FFT): ',num2str(kP_critical)])
% disp(['Real-frequency-bin(M-FFT): ',num2str(kvf)])
disp(['critically-sampled frequency upper-bound: ',num2str(fv_critical), ' Hz'])
disp(['critically-sampled frequency resolution: ',num2str(dfv_critical), ' Hz'])

%%
% x = zeros(1,L);
% for vv = 1:P
%     x = exp(1i*2*pi*kvf(vv)/M.*n) + x; 
% end
% x = exp(1i*2*pi*fv/srate.*n); 
x = exp(1i*2*pi*fv_sig/srate.*n);
disp(['Signal-frequency: ',num2str(fv_sig), ' Hz'])
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
% x_sysWin = x(1:M).* sys_win;
%
% figure(1),clf
% % subplot(211)
% hold on
% plot(0:M-1,real(x_sysWin),'b-','linew',2);
% % grid on
% title('Windowed-carrier-signal sequence''s realpart')
% % subplot(212)
% plot(0:M-1,imag(x_sysWin),'r-','linew',2);
% grid on
% 
% xlim([3000 5000])

%%
OSR = 2; % over-sampled ratio
chan_Nidx = kM/R+1; %% Caution: N channel decimation from M-FFT
dfv_osr = srate/M*(1 + 1/P*OSR);
disp(['Over-sampled frequency resolution: ',num2str(dfv_osr), ' Hz'])

Y_PFB = pfft(x,sys_win,N,R,OSR);
Y_samples = Y_PFB(chan_Nidx,:);
Z = fft(Y_samples(1:P));

disp('ends here')
%

%%
ksM = -M/2:M/2-1;
ksN = downsample(ksM,R);
ksP = -P/2:P/2-1; % FFT-bin index (normalized frequency) 

srate_PFB = srate/M;
frexM = ksM*srate/M;
frexN = downsample(frexM,R);
frexP = ksP*srate_PFB/P*OSR;

frexZ = kM*srate/M + frexP;

figure(3),clf
hold on
plot(frexN, fftshift(mag2db(abs(Y_PFB(:,1)))),'ko-','LineWidth',2)
grid on

figure(4),clf
hold on
plot(frexZ,fftshift(mag2db(abs(Z)./max(abs(Z)))),'bs-','LineWidth',2)
% ylim([-9 2])
% plot(ksP,(mag2db(abs(Z))),'bs','LineWidth',2)
% plot(mag2db(abs(Z)),'bs','LineWidth',2)
grid on
% ylim([-100 5])
%%
.