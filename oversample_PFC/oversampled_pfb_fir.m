function [fir_out] = oversampled_pfb_fir(x,h,M,R)
%oversampled_pfb_fir：多相滤波器的fir滤波器实现,2倍过采样率
% reg_out: 输出滤波器组输出的数据，[M,R*(W-1)+1]
% x: 输入信号
% h: 原型低通滤波器，长度为M*R
% M：多相滤波器的分支数 (后续FFT的点数，[0,2pi]划分的信道数)
% R：每一分支上滤波器的长度
%   

K = floor((length(x))/(M/2)); % 
% xp = reshape(x(1:M*R*W),[M W*R]);
x = x(1:M/2*K); % 截取数据长度为M/2*K
hp = reshape(h,[M R]); % 按矩阵M*R组织的LPF多相滤波器矩阵                                      

fir_out = zeros(M,K-R+1); % 存储每次load数据时，多相滤波器的输出
regp = reshape(x(1:M*R),[M R]); % 多相寄存器，初始值为x(1:M*R)"fully-loaded"
fir_out(:,1)=sum(hp.*regp,2); % PFB_fir第1列输出为fully-loaded的多相寄存器和多相滤波器的卷积

addr_ptr = M*R+1;
circ_buffer = zeros(1,M);

for cnt = 1:K-2*R % serpen_shift_load次数计数器 范围从1到K-2*R，对应x(MR+1:end)
    xs = x(addr_ptr:addr_ptr+M/2-1); %从输入数据中截取要载入的数据段。
    addr_ptr = addr_ptr+M/2; % 更新当前地址
    regp = serpen_shift_load(xs,regp);% 多相寄存器载入数据段
    fir_out(:,cnt+1)=sum(hp.*regp,2); % 计算PFB_fir相应列输出,从第2列开始，(列下标从1开始)
    
    %循环buffer操作,交换fir_out列的前M/2点和后M/2点，应用于fir_out的列2，4，6...
    if (mod(cnt+1,2)==0)
        circ_buffer(1:M/2)=fir_out(M/2+1:end,cnt+1);
        circ_buffer(M/2+1:end)=fir_out(1:M/2,cnt+1);
        fir_out(:,cnt+1)=circ_buffer;
    end
    
end
 
end

