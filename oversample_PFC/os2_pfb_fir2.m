function [fir_out] = os2_pfb_fir2(x,h,M,R)
%oversampled_pfb_fir：多相滤波器的fir滤波器实现,2倍过采样率
% fir_out: 输出滤波器组输出的数据，[M,(K-2R)]
% x: 输入信号
% h: 原型低通滤波器，长度为M*R,
% M：多相滤波器的分支数 (后续FFT的点数，[0,2pi]划分的信道数)
% R：每一分支上滤波器的长度
%   

K = floor((length(x))/(M/2)); % 
% xp = reshape(x(1:M*R*W),[M W*R]);
x_1st = x(1:M/2*(K-1)); % 原序列截取数据长度为M/2*(K-1)
x_2nd = x(M/2+1:M/2*K); % 原序列移位M/2点后截取数据长为M/2*(K-1)

hp = reshape(h,[M R]); % 按矩阵M*R组织的LPF多相滤波器矩阵                                      
% hp = (upsample(hp',2))'; % 在多相滤波器矩阵的每一分支上插入寄存器

fir_out = zeros(M,K-1); % 存储每次load数据时，多相滤波器的输出
regp = zeros(M,2*R-1);
% regp = cat(1,reshape(x_2nd(1:M*R),[M/2 2*R]),reshape(x_1st(1:M*R),[M/2 2*R])); % 多相寄存器，初始值为"fully-loaded"
% regp_d = (downsample(regp',2))';
% fir_out(:,1)=sum(hp.*regp_d,2); % PFB_fir第1列输出为fully-loaded的多相寄存器和多相滤波器的卷积

fp=1;
circ_buffer = zeros(M,1);

for cnt = 1:K-1 % load了K-1次,计数器范围从1到K-1
    xs1 = flipud(reshape(x_1st(fp:fp+M/2-1),[M/2,1])); %从输入数据1中截取要载入的数据段。
    xs2 = flipud(reshape(x_2nd(fp:fp+M/2-1),[M/2,1])); %从输入数据2中截取要载入的数据段。
    xs = cat(1,xs2, xs1);
    fp = fp+M/2; % 更新当前地址
    regp(:,2:end) = regp(:,1:end-1);% 多相寄存器从左到右载入数据段
    regp(:,1)=xs;
    regp_1 = (downsample((fliplr(regp))',2))';
    regp_2=fliplr(regp_1);
    fir_out(:,cnt)=sum(hp.*regp_2,2); % 计算PFB_fir相应列输出,从第1列开始，(列下标从1开始)
    
    %循环buffer操作,交换fir_out列的前M/2点和后M/2点，应用于fir_out的列2，4，6...
    if (mod(cnt,2)==1)
        circ_buffer(1:M/2)=fir_out(M/2+1:end,cnt);
        circ_buffer(M/2+1:end)=fir_out(1:M/2,cnt);
        fir_out(:,cnt)=circ_buffer;
    end
    
end
 
end

