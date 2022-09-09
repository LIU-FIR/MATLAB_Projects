function [Rp_bin] = pfftBinResponse(win,N,R,dks,kvM,kvspan)
% polyphase_FFT输出某个bin上的响应
% Rp: poly phase FFT outputs(N*Frames)，N-channel in each column,
%               Frames depends on OSR(Over-sampled-ratio)
% dks: k bins fractional resolution.
% kvM: selected channel from M-channel of pfft output, M = N*R
% R: pfft taps;
M = N*R;
kvN = kvM/R;
chan_Nidx = kvN + 1;
n = 0:M-1;
ks = kvM-kvspan:dks:kvM+kvspan; %
% disp(num2str(ks))

Rp_bin = zeros(length(ks),1);
for kk = 1:length(ks)

    x_kk = exp(1i*2*pi*ks(kk)./M.*n); % generate a signal with new freq.
    
    % compute the pfft frames.
    Y_PFB = pfft(x_kk,win,N,R,1);
    
    Rp_bin(kk,:)=abs(Y_PFB(chan_Nidx,:)); % pick up the kv-th bin in pfft's bin. 
end

end

