function reg_out = serpen_shift_load(vcol,reg_mat)
%UNTITLED 此处显示有关此函数的摘要
%   vcol: M/2*1 向量，来自于要载入的数据
%   reg_mat: M*R 矩阵，不断从外部数据载入vcol，vcol载入口为右下角，
[M,R] = size(reg_mat);
reg_out = zeros(M,R);
reg_out(1:M/2,:)=reg_mat(M/2+1:end,:);
reg_out(M/2+1:end,1:end-1)=reg_mat(1:M/2,2:end);
reg_out(M/2+1:end,end)=vcol;
end

