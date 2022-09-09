%%
% 验证长序列 M-pnt-PFFT的streamlining spectral response
clear;clc;
L = 20000; % original signal length.
M = 256; % polyphase-FFT-length
R = 8; % Taps in each polyphase branch. 
N = M*R; % Polyphase Signal sequence block size
n = 0:L-1;   % sampling index (time)
k = -N/2:N/2-1; % FFT-bin index (normalized frequency) 
srate = 1e3; % sampling rate.

ham_win = (hamming(N))'; % hamming window function
sinc_win = fir1(N-1,1/M,rectwin(N)); % sinc window function
sinc_win_norm = sinc_win./max(sinc_win); % nomralized sinc window function
sys_win = sinc_win_norm.*ham_win; % composed system window function

%%
os_ratio_str = '4/3';
os_ratio = str2num(os_ratio_str); % the oversampled ratio. os_ratio>1 means that outputs from PFFT's bins refresh faster.
slide_len = N/os_ratio; % when being loaded,the start position in orig.sequence increase by slide_len

load_times = floor((L-N)/slide_len)+1; % deducted from geometric relation
disp(['load_times: ' num2str(load_times)])

Yp_frames = zeros(M,load_times);

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
% subplot(121)
hold on

for ip = 1:load_times
    plot(ks, mag2db(Rp_frex(:,ip)/max(Rp_frex(:,ip)))+1.5*(ip-1),'LineWidth',1)
end
set(gca,'ytick',[])
xlabel('PFFT bins')
ylabel('PFFT power (dB)')
title(['Bin power response around ', num2str(kv),'-th PFFT-bin of ',num2str(load_times),' frames'])

%%
figure(4),clf
% subplot(211)
% imagesc(abs(Yp_frames))
% axis xy
% % axis square
% colormap('hot')
% colorbar
% xlabel('PFFT frames (in time)')
% ylabel('PFFT bins')
% title(['Output streams of PFFT; the signal at ',num2str(kv),'-th bin'])

% subplot(212)
imagesc(Rp_frex)
axis xy
% axis square
colormap('hot')
colorbar
set(gca,'ytick',[101 201 301 401 501],'yticklabel',{'14','15','16','17','18'})
xlabel('PFFT time frames (m)')
ylabel('PFFT bin-freq.')
title(['Bin power-response around ', num2str(kv),'-th PFFT-bin of ',num2str(load_times),' frames'])
%% done
