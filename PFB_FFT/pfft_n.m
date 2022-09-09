function [y] = pfft_n(x,win,N,R,dsi)
% polyphase_fft的多相结构实现
% pfft:poly phase FFT 输出的数据，[M,R*(W-1)+1]
% x: 输入信号
% win: 窗函数，长度为M*R
% M：多相FFT的通道数 
% R：每一分支上窗函数的长度
% dsi: M*R点FFT的抽头位置(down-sampled-index)，默认为0，范围是：[0 R-1]

W= floor(length(x)/N);
winp = reshape(win,[N R]); % 按矩阵M*R组织的LPF系数

reg = zeros(1,N*R);% PFFT寄存器，M*R-length

K = W-R+1;
y = zeros(N,K); % 存储每次load数据时，多相FFT的输出

N = N*R;
m = (0:N-1)';
pha_m = exp(-1i*2*pi*m*dsi/N); % 抽头位置应用在通道m上的相位改正

r = 0:R-1;
pha_r = exp(-1i*2*pi*r*dsi/R); % 抽头位置应用在所有通道r-th tap上的相位改正
winp = winp.*pha_r;

for cnt = 1:K %fp: frame pointer, 范围从1到R*W
    if cnt==1
        reg = x(1:N*R);
        regp = reshape(reg,[N R]);
        y(:,1)=fft(sum(winp.*regp,2).*pha_m);
    end
    if K==1
        break;
    end
    reg(1:N*R-N)=reg(N+1:end);% 多相寄存器移位载入数据段
    reg(N+1:end)= x(cnt*N+1:cnt*N+N*(R-1));%从输入数据流中截取载入的数据段
    regp = reshape(reg,[N R]);
    y(:,cnt)=fft(sum(winp.*regp,2).*pha_m);
 
end

end

