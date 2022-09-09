%%
% 验证长序列FFT多相结构的连续输出
clear;clc;
L = 20000; % original signal length.
M = 256; % polyphase-FFT-length
R = 8; % Taps in each polyphase branch. 
N = M*R; % Polyphase Signal sequence block size
n = 0:L-1;   % sampling index (time)
k = -N/2:N/2-1; % FFT-bin index (normalized frequency) 
srate = 1e3; % sampling rate.
f1 = 100; % signal frequency.
f2 = 200;
% kv1 = 16; % signal frequency bin in M-pnt FFT.
% kv2 = 32; % signal frequency bin in M-pnt FFT.
kv = 32; % signal frequency integer-part for M-pnt FFT.
dkv = 0.05; % signal frequency fraction-part for M-pnt FFT.
kvf = kv + dkv;
% complex sinsoid signal (input)
% x = exp(1i*2*pi*f1/srate.*n)+1.5*exp(1i*2*pi*f2/srate.*n); 
% pha1 = pi/16;
% pha2 = 3*pi/16;
% x = exp(1i*2*pi*kv1/M.*n + 1i*pha1)+1.5*exp(1i*2*pi*kv2/M.*n + 1i*pha2); 
pha = 0;
x = exp(1i*2*pi*kvf/M.*n + 1i*pha); 
% win = (rectwin(N))'; 
ham_win = (hamming(N))'; % hamming window function
sinc_win = fir1(N-1,1/M,rectwin(N)); % sinc window function
sinc_win_norm = sinc_win./max(sinc_win); % nomralized sinc window function
sys_win = sinc_win_norm.*ham_win; % composed system window function

x_win = sys_win.*x(1:N);
X_win = fft(x_win);
X_win_mag = abs(X_win);

% plot time squence data.
figure(1),clf
subplot(2,1,1)
plot(0:N-1,real(x(1:N)),'k-');
ylim([0 3])
xlim([500 1500])
grid on
title('Orig. signal sequence')

subplot(2,1,2)
plot(0:N-1,abs(x_win),'b-');
grid on
title('Windowed signal sequence')
ylim([0 3])
xlim([500 1500])

figure(2),clf
plot(k,fftshift(X_win_mag),'k-');
grid on
title('Windowed Signal spectrum')
xlim([0 500])

%%
os_ratio = 2; % the oversampled ratio. os_ratio>1 means that outputs from PFFT's bins refresh faster.
slide_len = N/os_ratio; % when being loaded,the start position in orig.sequence increase by slide_len

load_times = floor((L-N)/slide_len)+1; % deducted from geometric relation
disp(['load_times: ' num2str(load_times)])

Yp_frames = zeros(M,load_times);
for cnt = 1:load_times
    % for every sliding, compute the current head and tail of the
    % ready-to-pfft sequence.
    head_ptr = 1+(cnt-1)*slide_len;
    tail_ptr = 1+(cnt-1)*slide_len + N-1;
    
    x_reg = x(head_ptr:tail_ptr); % load the ready-to-pfft sequence in register.
    Yp_frames(:,cnt) = pfft(x_reg,sys_win,M,R);
end


Yp_mag_frames = abs(Yp_frames);
% Yp_kv1_frames = Yp_frames(kv1+1,:);
% Yp_kv2_frames = Yp_frames(kv2+1,:);
% angles_kv1 = angle(Yp_kv1_frames);
% angles_kv2 = angle(Yp_kv2_frames);

Yp_kv_frames = Yp_frames(kv+1,:);

m = 0:load_times-1;
phasor_lags = pi*m*dkv*R;
angles_kv = unwrap(angle(Yp_kv_frames));

% Yp_kv_frames_rev = Yp_kv_frames.*(-1).^(-m*kv).*exp(-1i*pi*m*dkv);
angles_kv_rev = angles_kv-phasor_lags;

figure(3),clf
hold on
plot(m,angles_kv,'ko-','linew',1,'markerfacecolor','b')
plot(m,angles_kv_rev,'ro-','linew',1,'markerfacecolor','r')
% plot(m,phasor_lags,'b-','linew',1)
%%
% plot all PFFT frames.1) at kv1-bin and kv2-bin; 2) from all fft-bins
figure(2),clf
subplot(121)
hold on
plot(Yp_mag_frames(kv1+1,:),'k-','linew',2)
plot(Yp_mag_frames(kv2+1,:),'b-','linew',2)
% plot(Yp_mag_frames(120,:),'m-','linew',2)
legend({[num2str(kv1) '-th bin'] [num2str(kv2) '-th bin']},'location','southeast')
yrange = get(gca,'ylim');
set(gca,'ylim',[0 yrange(2)])
subplot(122)
imagesc(Yp_mag_frames)
axis xy
colormap('jet')
colorbar
%%
% Compute spectral response from a given pfft-bin
kv = 16;
ks = kv-3:0.01:kv+3; % frequencies to be computed.
% tN = (0:N-1); % original N-pnt length sample indices
tL = (0:L-1); % original L-pnt length sample indices
% tM = (0:M-1); % M-pnt length sample indices.

% Rp_frex: polyphase-FFT freq.response from the N-pnt sequence.  
Rp_frex = zeros(length(ks),load_times);
for kk = 1:length(ks)

    xk = exp(1i*2*pi*ks(kk)./M.*tL); % generate a signal with new freq.
    
    % compute the pfft frames.
    for cnt = 1:load_times
        % for every sliding, compute the current head and tail of the
        % ready-to-pfft sequence.
        head_ptr = 1+(cnt-1)*slide_len;
        tail_ptr = 1+(cnt-1)*slide_len + N-1;

        xk_reg = xk(head_ptr:tail_ptr); % load the ready-to-pfft sequence in register.
        Yp_frames(:,cnt) = pfft(xk_reg,sys_win,M,R);
    end    
    
    Rp_frex(kk,:)=abs(Yp_frames(kv+1,:)); % pick up the kv-th bin. 

end
disp('ends here')
%%
figure(3),clf
hold on
for ip = 1:load_times
    plot(ks, mag2db(Rp_frex(:,ip)/max(Rp_frex(:,ip)))+(ip-1),'LineWidth',1)
end


%%
ham_win_pmat = reshape(ham_win,[M R]);% hamming window polyphase matrix
sys_win_pmat = reshape(sys_win,[M R]);% system window polyphase matrix
x_pmat = reshape(x,[M R]); % signal polyphase matrix
z = sum(ham_win_pmat.*(x_pmat),2);
Yp0 = fft(z); % use matrix-method to compute polyphase FFT.
Yp= pfft(x,ham_win,M,R); % use pfft() to compute polyphase FFT.
X_ds = downsample(X_win,R);% resample the orignal long FFT.

kd = (0:M-1);


figure(2),clf
plot(kd,abs(Yp0),'ko')
hold on
plot(kd,abs(Yp),'r*')
hold on
plot(kd,abs(X_ds),'b^')
grid on
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