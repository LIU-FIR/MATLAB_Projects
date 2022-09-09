%%
% 验证coarse M-pnt-pfft_cs + fine P-pnt-FFT
clear;clc;

N = 2048; % polyphase-FFT-length
R = 4; % Taps in each polyphase branch. 
M = N*R; % Polyphase Signal sequence block size
P = 64; % 2-nd fine FFT-lingth
L = M*120; % original signal length.
n = 0:L-1;   % sampling index (time)
kM = -M/2:M/2-1; % FFT-bin index (normalized frequency) 
srate = 2048e6; % sampling rate.

os_ratio_str = '2';
os_ratio = str2num(os_ratio_str); % the oversampled ratio. os_ratio>1 means that outputs from pfft_cs's bins refresh faster.
slide_len = M/os_ratio; % when being loaded,the start position in orig.sequence increase by slide_len
num_frames = floor((L-M)/slide_len)+1; % deducted from geometric relation
disp(['num_frames: ' num2str(num_frames)])

kvM = 128; % signal frequency integer-part for M-pnt pfft_cs.

% set 2-nd FFT bin,equivalent to fraction-part for M-pnt pfft_cs.
rvs = -P/2:P/2-1;
% rv1 = 1;
% rv2 = 7;
% rv3 = 12;
% dkv1 = rv1/P/R * os_ratio; 
% dkv2 = rv2/P/R * os_ratio; 
% dkv3 = rv3/P/R * os_ratio; 
% x = exp(1i*2*pi*k/M.*n + 1i*pha)...
%     + exp(1i*2*pi*kvf1/M.*n + 1i*pha)...
%     + exp(1i*2*pi*kvf2/M.*n)...
%     + exp(1i*2*pi*kvf3/M.*n); 

dks = rvs/P * os_ratio;
% siganal frequency (integar and fraction) for M-pnt pfft_cs.
kvs = kvM + dks;
pha = pi/16;

x = zeros(1,L);
for cc = 1:P
    x = x + exp(1i*2*pi*kvs(cc)/M.*n);
end
%%
rect_win = (rectwin(M))';
ham_win = (hamming(M))'; % hamming window function
sinc_win = fir1(M-1,1/N,rectwin(M)); % sinc window function
% sinc_win_norm = sinc_win./max(sinc_win); % nomralized sinc window function
sinc_win_norm = sinc_win./sum(sinc_win); % nomralized sinc window function


% sinc_win1 = sinc(pi.*((0:M-1)-M/2)./(4*N))/4;
% sinc_win1 = sinc_win1./sum(sinc_win1);
% 
% figure(1),clf
% hold on
% plot(sinc_win,'k-','linew',1)
% plot(sinc_win1,'b-','linew',1)
% 
% y_sincWin = fft(sinc_win);
% y1_sincWin = fft(sinc_win1);
% 
% figure(2),clf
% hold on 
% plot(fftshift(mag2db(abs(y_sincWin))),'k-','linew',1)
% plot(fftshift(mag2db(abs(y1_sincWin))),'b-','linew',1)

%%
sys_win = sinc_win_norm.*ham_win; % composed system window function
% sys_win = rect_win;

x_win = sys_win.*x(1:M);
X_win = fft(x_win);
X_win_mag = abs(X_win);

% plot time squence data.
figure(1),clf
subplot(2,1,1)
plot(0:M-1,real(x(1:M)),'k-','linew',1);
% ylim([0 1.5])
xlim([500 1500])
grid on
title('Orig. signal sequence')

subplot(2,1,2)
plot(0:M-1,abs(x_win),'b-','linew',1);
grid on
title('Windowed signal sequence')
% ylim([0 1.5])
xlim([500 1500])

%%

Yp_frames = zeros(N,num_frames);
Yp_rect_frames = zeros(N,num_frames);
for cnt = 1:num_frames
    % for every sliding, compute the current head and tail of the
    % ready-to-pfft_cs sequence.
    head_ptr = 1+(cnt-1)*slide_len;
    tail_ptr = 1+(cnt-1)*slide_len + M-1;
    
    x_reg = x(head_ptr:tail_ptr); % load the ready-to-pfft_cs sequence in register.
    Yp_frames(:,cnt) = pfft_cs(x_reg,sys_win,N,R);
    Yp_rect_frames(:,cnt) = pfft_cs(x_reg,rect_win,N,R);
end


Yp_kvN_frames = Yp_frames(kvM/R+1,1:P); % kvM->kvN, N-chan decimated from M-FFT
Yp_kvN_rect_frames = Yp_rect_frames(kvM/R+1, 1:P);

Z = fft(Yp_kvN_frames);
Z_rect = fft(Yp_kvN_rect_frames);
r = -P/2:P/2-1;

Z_db = mag2db(abs(Z)./abs(Z(1)));
Z_rect_db = mag2db(abs(Z_rect)./abs(Z_rect(1)));
disp('ends here')
%%

figure(2),clf
hold on
% stem(r,mag2db(fftshift(abs(Z))),'ks','markerfacecolor','b')
plot(r,fftshift(Z_db),'bo-','markerfacecolor','b','linew',2)
plot(r,fftshift(Z_rect_db),'ks-','markerfacecolor','k','linew',2)
grid on
xlim([-P/2-2 P/2+2])
ylim([-9 2])
% patch('xdata',[-16 -16 15 15],'ydata',[-9 0 0 -9],'edgecolor','m','facecolor','none','linew',1)

legend({'hanning*sinc-window','rect-window'})
xlabel('FFT bins')
ylabel('Power response: (dB)')
title(['Fine ',num2str(P),'pnt-FFT bin response: ','os\_rate= ',os_ratio_str,'; rect-win v.s. hanning*sinc-win'])
% title(['2-nd ',num2str(P),'pnt-FFT spectral response: ','oversampled with ',os_ratio_str,'; hanning*sinc-win'])
%%
% plot(0,mag2db(abs(Z(1))./abs(Z(1))),'ro','markerfacecolor','r')
% plot(rv1,mag2db(abs(Z(rv1+1))./abs(Z(1))),'ro','markerfacecolor','r')
% plot(rv2,mag2db(abs(Z(rv2+1))./abs(Z(1))),'ro','markerfacecolor','r')
% plot(rv3,mag2db(abs(Z(rv3+1))./abs(Z(1))),'ro','markerfacecolor','r')
% % xticks([-P/2 0 rv1 rv2 rv3 P/2-1])