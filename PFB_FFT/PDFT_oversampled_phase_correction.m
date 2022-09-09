%%
% L-pnt x(n) sequence:
% M-PFFT oversampling_rate = 2, 4/3
% phase consistent correction.

clear;clc;
L = 20000; % original signal length.
M = 256; % polyphase-FFT-length
R = 8; % Taps in each polyphase branch. 
N = M*R; % Polyphase Signal sequence block size
n = 0:L-1;   % sampling index (time)
k = -N/2:N/2-1; % FFT-bin index (normalized frequency) 
srate = 1e3; % sampling rate.
kv = 32; % signal frequency integer-part for M-pnt FFT.
dkv = 0.05; % signal frequency fraction-part for M-pnt FFT.
% dkv = 0;
kvf = kv + dkv;

pha = pi/16;
x = exp(1i*2*pi*kvf/M.*n + 1i*pha); 

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
os_ratio_str = '1/1'; % the oversampled ratio. os_ratio>1 means that outputs from PFFT's bins refresh faster.
os_ratio = str2num(os_ratio_str);
slide_len = N/os_ratio; % sliding length to form a new frame.

% number of consecutive frames from the orig.squence, 
% deducted from geometric relation.
num_frames = floor((L-N)/slide_len)+1; 
disp(['num_frames: ' num2str(num_frames)])
disp(['slide_length: ' num2str(slide_len)])
disp(['os_ratio: ' num2str(os_ratio)])
%
Yp_frames = zeros(M,num_frames);

for cnt = 1:num_frames
    % for every sliding, compute the current head and tail of the
    % ready-to-pfft sequence.
    head_ptr = 1+(cnt-1)*slide_len;
    tail_ptr = 1+(cnt-1)*slide_len + N-1;
    
    x_reg = x(head_ptr:tail_ptr); % load the ready-to-pfft sequence in register.
    Yp_frames(:,cnt) = pfft(x_reg,sys_win,M,R);
end


% Yp_mag_frames = abs(Yp_frames);

Yp_kv_frames = Yp_frames(kv+1,:);

m = 0:num_frames-1;

phase_rev = -2*pi*R*dkv/os_ratio.*m;
angles_kv = unwrap(angle(Yp_kv_frames));

% Yp_kv_frames_rev = Yp_kv_frames.*(-1).^(-m*kv).*exp(-1i*pi*m*dkv);
angles_kv_rev = angles_kv + phase_rev;

figure(3),clf
hold on
title(['Phase variations in M-PFFT streamline: ','os\_rate = ',os_ratio_str,...
    '; R\cdot\delta=',num2str(R*dkv),'; N=',num2str(N),'; M=',num2str(M)])
plot(m,angles_kv,'ks-','linew',1,'markerfacecolor','b','markersize',10)
plot(m,angles_kv_rev,'ro-','linew',1,'markerfacecolor','r','markersize',10)
xlabel('Sampling points at the k-th bin.(a.u.)')
ylabel('Phase angle (rad)')
legend({'Before correction';'After correction'},'location','northwest')
grid on
% ylim([0 2])
% plot(m,phasor_lags,'b-','linew',1)
%%


%% done
