%%
% 验证长序列FFT多相结构的连续输出
clear;clc;

M = 8192; % Long-FFT-length
R = 4; % Taps in each polyphase branch. 
N = M/R; % Polyphase FFT-length
L = M; % original signal length.
n = 0:L-1;   % sampling index (time)
kM = -M/2:M/2-1; % FFT-bin index (normalized frequency) 
srate = 2048e6; % sampling rate.

kv = 160; % signal frequency integer-part for M-pnt long-FFT.
% dkv = 0.15; % signal frequency fraction-part for M-pnt FFT.
dkv = 0;
kvf = kv + dkv;

x = exp(1i*2*pi*kvf/M.*n); 
% win = (rectwin(N))'; 
ham_win = (hamming(M))'; % hamming window function
ham_win = ham_win./sum(ham_win);
% sinc_win = fir1(N-1,1/M,rectwin(N)); % sinc window function

% sinc_win = sinc(pi.*((0:M-1)-M/2)./N);
% sinc_win = sinc(pi.*((0:M-1)-M/2)./(2*N))/2;
sinc_win = sinc(pi.*((0:M-1)-M/2)./(4*N))/4;
sinc_win = sinc_win./sum(sinc_win);
%%
sys_win = sinc_win.*ham_win; % composed system window function
sys_win = sys_win./sum(sys_win);
rect_win = ones(1,M)./M;

x_sysWin = sys_win.*x(1:M);
X_sysWin = fft(x_sysWin)./M;
X_sysWinMag = abs(X_sysWin);

X_sincWin = fft(sinc_win.*x(1:M))./M;
X_sincWinMag = abs(X_sincWin);

X_hamWin = fft(ham_win.*x(1:M))./M;
X_hamWinMag = abs(X_hamWin);

X_recWin = fft(rect_win.*x(1:M))./M;
X_recWinMag = abs(X_recWin);
%%
% plot time squence data.
figure(1),clf
% subplot(211)
hold on
plot(0:M-1,real(x_sysWin),'b-','linew',2);
% grid on
title('Windowed-carrier-signal sequence''s realpart')
% subplot(212)
plot(0:M-1,imag(x_sysWin),'r-','linew',2);
grid on
% title('Windowed-carrier-signal sequence''s imagpart')
% title('Windowed signal sequence''s imagpart')

% ylim([0 3])
xlim([3000 5000])

figure(2),clf
hold on
plot(kM,fftshift(mag2db(X_recWinMag./max(X_recWinMag))),'r-','linew',2);
plot(kM,fftshift(mag2db(X_sincWinMag./max(X_recWinMag))),'k-','linew',2);
plot(kM,fftshift(mag2db(X_hamWinMag./max(X_recWinMag))),'m-','linew',2);
plot(kM,fftshift(mag2db(X_sysWinMag./max(X_recWinMag))),'b-','linew',2);
grid on
title('Windowed-carrier-signal spectrum')
legend({'Rect-wined signal-FFT mag','Sinc-wined signal-FFT mag',...
    'Ham-wined signal-FFT mag','Sys-wined signal-FFT mag'});
xlim([0 300])
xlabel('FFT bins')
ylabel('Power (dB)')

%%
os_ratio = 1; % the oversampled ratio. os_ratio>1 means that outputs from PFFT's bins refresh faster.
slide_len = M/os_ratio; % when being loaded,the start position in orig.sequence increase by slide_len

load_times = floor((L-M)/slide_len)+1; % deducted from geometric relation
disp(['load_times: ' num2str(load_times)])

% plot(m,phasor_lags,'b-','linew',1)
%
% Compute M-long-FFT over all frames.
% X_frames = zeros(M,load_times);
% ks = kv-10:0.1:kv+10; % frequencies to be computed.
% Rp_frex_MFFT = zeros(length(ks),load_times);
% for kk = 1:length(ks)
% 
%     x_kk = exp(1i*2*pi*ks(kk)./M.*n); % generate a signal with new freq.
%     x_frames = reshape(x_kk(1:M*load_times),[M,load_times]);    
%     X_frames = fft(x_frames.*sys_win');
%     Rp_frex_MFFT(kk,:)=abs(X_frames(kv+1,:)); % pick up the kv-th bin in pfft's bin. 
% 
% end
% disp('computing ends')
% %
% figure(3),clf
% hold on
% plot(ks, mag2db(Rp_frex_MFFT(:,1)./max(Rp_frex_MFFT(:,1))),'LineWidth',1)
% ylim([-100 0])


%%
% plot all PFFT frames.1) at kv1-bin and kv2-bin; 2) from all fft-bins
%
% Compute spectral response from a given pfft-bin
kvspan = 10;
Yp_rectWin_frames = zeros(N,load_times);
Yp_hamWin_frames = zeros(N,load_times);
Yp_sincWin_frames = zeros(N,load_times);
Yp_sysWin_frames = zeros(N,load_times);
ks = kv-kvspan:0.01:kv+kvspan; % frequencies to be computed.

% ks = kv;
% tL = (0:L-1); % original L-pnt length sample indices
% tM = (0:M-1); % M-pnt length sample indices.

% Rp_frex: polyphase-FFT freq.response from the N-pnt sequence.  
Rp_rectWin_NPFFT = zeros(length(ks),load_times);
Rp_hamWin_NPFFT = zeros(length(ks),load_times);
Rp_sincWin_NPFFT = zeros(length(ks),load_times);
Rp_sysWin_NPFFT = zeros(length(ks),load_times);
RpL_sysWin_NPFFT = zeros(length(ks),load_times);
RpR_sysWin_NPFFT = zeros(length(ks),load_times);
for kk = 1:length(ks)

    x_kk = exp(1i*2*pi*ks(kk)./M.*n); % generate a signal with new freq.
    
    % compute the pfft frames.
    for cnt = 1:load_times
        % for every sliding, compute the current head and tail of the
        % ready-to-pfft sequence.
        head_ptr = 1+(cnt-1)*slide_len;
        tail_ptr = 1+(cnt-1)*slide_len + M-1;

        xk_reg = x_kk(head_ptr:tail_ptr); % load the ready-to-pfft sequence in register.
        Yp_rectWin_frames(:,cnt) = pfft_cs(xk_reg,rect_win,N,R);
%         Yp_hamWin_frames(:,cnt) = pfft(xk_reg,ham_win,N,P);
%         Yp_sincWin_frames(:,cnt) = pfft(xk_reg,sinc_win,N,P);
        Yp_sysWin_frames(:,cnt) = pfft_cs(xk_reg,sys_win,N,R);
    end    
    
    Rp_rectWin_NPFFT(kk,:)=abs(Yp_rectWin_frames(kv/R+1,:)); % pick up the kv-th bin in pfft's bin. 
%     Rp_hamWin_NPFFT(kk,:)=abs(Yp_hamWin_frames(kv/P+1,:)); % pick up the kv-th bin in pfft's bin. 
%     Rp_sincWin_NPFFT(kk,:)=abs(Yp_sincWin_frames(kv/P+1,:)); % pick up the kv-th bin in pfft's bin. 
    Rp_sysWin_NPFFT(kk,:)=abs(Yp_sysWin_frames(kv/R+1,:)); % pick up the kv-th bin in pfft's bin. 
    RpR_sysWin_NPFFT(kk,:)=abs(Yp_sysWin_frames(kv/R+1+1,:)); % pick up the kv-th bin in pfft's bin. 
    RpL_sysWin_NPFFT(kk,:)=abs(Yp_sysWin_frames(kv/R+1-1,:)); % pick up the kv-th bin in pfft's bin. 

end
disp('ends here')

%%
figure(4),clf
hold on
plot(ks, mag2db(Rp_rectWin_NPFFT(:,1)./max(Rp_rectWin_NPFFT(:,1))),'c-','LineWidth',2)
% plot(ks, mag2db(Rp_hamWin_NPFFT(:,2)./max(Rp_rectWin_NPFFT(:,2))),'m-','LineWidth',2)
% plot(ks, mag2db(Rp_sincWin_NPFFT(:,2)./max(Rp_rectWin_NPFFT(:,2))),'k-','LineWidth',2)
plot(ks, mag2db(Rp_sysWin_NPFFT(:,1)./max(Rp_rectWin_NPFFT(:,1))),'b-','LineWidth',2)
ylim([-100 5])
%%
figure(5),clf
hold on
plot(ks, mag2db(Rp_sysWin_NPFFT(:,1)./max(Rp_rectWin_NPFFT(:,1))),'b-','LineWidth',2)
% plot(ks, mag2db(Rp_hamWin_NPFFT(:,2)./max(Rp_rectWin_NPFFT(:,2))),'m-','LineWidth',2)
% plot(ks, mag2db(Rp_sincWin_NPFFT(:,2)./max(Rp_rectWin_NPFFT(:,2))),'k-','LineWidth',2)
plot(ks, mag2db(RpR_sysWin_NPFFT(:,1)./max(Rp_rectWin_NPFFT(:,1))),'r-','LineWidth',2)
plot(ks, mag2db(RpL_sysWin_NPFFT(:,1)./max(Rp_rectWin_NPFFT(:,1))),'k-','LineWidth',2)
ylim([-100 5])
grid on
%%
Rp1_sysWin = pfftBinResponse(sys_win,N,R,0.01,kv,kvspan);

figure(6)
hold on
plot(ks, mag2db(Rp_sysWin_NPFFT(:,1)./max(Rp_rectWin_NPFFT(:,1))),'bo','LineWidth',2)
plot(ks, mag2db(Rp1_sysWin./max(Rp1_sysWin)),'k.','LineWidth',2)