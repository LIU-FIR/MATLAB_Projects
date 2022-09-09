function [circ_buffer] = oversampled_pfb_fir1(x,h,M,R)
%oversampled_pfb_fir：多相滤波器的fir滤波器实现,2倍过采样率
% reg_out: 输出滤波器组输出的数据，[M,R*(W-1)+1]
% x: 输入信号
% h: 原型低通滤波器，长度为M*R
% M：多相滤波器的分支数 (后续FFT的点数，[0,2pi]划分的信道数)
% R：每一分支上滤波器的长度
%   
Ri = 2*R; % 每一分支上内插后滤波器的长度
W = floor(length(x)/(M/2*Ri)); % 按M/2通道排列和分支多相滤波器长度为R组织，计算信号“块”数。
x = x(1:M/2*Ri*W);
xp_half = reshape(x,[M/2 W*Ri]);
xp = cat(1,xp_half,xp_half);

hp = reshape(h,[M R]); % 按矩阵M*R组织的LPF系数
hp_i = zeros(M,Ri); % 内插后的LPF系数

%对每一分支的滤波器系数进行2倍内插
for j = 1:M/2
    hp_i(j,:) = upsample(hp(j,:),2);
    hp_i(j+M/2,:) = upsample(hp(j,:),2,1);
end
% regp = zeros(M,R);
ind_m = Ri*(W-1)+1; % load最后一块数据(M*R)时，数据的帧索引
reg_out = zeros(M,ind_m); % 存储每次load数据时，多相滤波器的输出
for fp = 1:ind_m %fp: frame pointer, 范围从1到ind_m
    regp = xp(:,fp:fp+Ri-1);
    reg_out(:,fp)=sum(hp_i.*regp,2);
end
 
%循环buffer,交换前M/2和后M/2的通道
circ_buffer = zeros(M,ind_m);
circ_buffer(1:M/2,:)=reg_out(M/2+1:end,:);
circ_buffer(M/2+1:end,:)=reg_out(1:M/2,:);
% circ_buffer = reg_out;
end

