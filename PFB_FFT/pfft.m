function [Y_PFB_frames] = pfft(x,win,N,R,OSR)
% polyphase_fft: M=N*R, N-channel decimated from M-FFT
% Y_PFB_frames: poly phase FFT outputs(N*Frames)，N-channel in each column,
%               Frames depends on OSR(Over-sampled-ratio)
% x: L-pnt input signal
% win: N*R window function
% M：pfft channels 
% R：pfft taps
% OSR: Over-sampled-ratio

L = length(x); % signal's length, L > N*R
M = N*R; % M-long-FFT points
winp = reshape(win,[N R]); % 

slide_len = M/OSR; % when being loaded,the start position in orig.sequence steps forward by slide_len

load_times = floor((L-M)/slide_len)+1; % deducted from geometric relation

Y_PFB_frames = zeros(N,load_times); % Poly-phase-FFT-filterBank outputs, N-channel in each column.

for cnt = 1:load_times
    % for every sliding, compute the current head and tail of the
    % ready-to-pfft_cs sequence.
    head_ptr = 1+(cnt-1)*slide_len;
    tail_ptr = 1+(cnt-1)*slide_len + M-1;

    xk_reg = x(head_ptr:tail_ptr); % load the ready-to-pfft_cs sequence in register.
    regp = reshape(xk_reg,[N R]);
 
    Y_PFB_frames(:,cnt) = fft(sum(winp.*regp,2));
end    

end

