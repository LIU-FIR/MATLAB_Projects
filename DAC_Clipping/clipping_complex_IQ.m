%%
clear;clc;

%% 
% 形成IQ两路的载波信号
N = 1024;
n = 0:N-1;
fs = 320e6; % 采样率

ks = [20 64 100 120]; % 载波对应的bins
phases = randperm(N,4)/N*2*pi;

Nc = 4;

AF = 2; % DAC full-scale
P_avg = AF.^2; %I/Q的复信号总功率
A = AF/Nc^0.5; 

%%
ind_ks = ks+1;

Xs = zeros(N,Nc);

for vv = 1:Nc
    Xs(ks(vv)+1,vv) = A*exp(1i*phases(vv));
end

s_arr = N*ifft(Xs);
s_sum = sum(s_arr,2);
X_sum = 1/N*fft(s_sum);

s_DAC_clip_real = real(s_sum);
s_DAC_clip_imag = imag(s_sum);
s_DAC_clip_real(s_DAC_clip_real>AF) = AF;
s_DAC_clip_real(s_DAC_clip_real<-AF) = -AF;
s_DAC_clip_imag(s_DAC_clip_imag>AF) = AF;
s_DAC_clip_imag(s_DAC_clip_imag<-AF) = -AF;
s_DAC_clip = complex(s_DAC_clip_real,s_DAC_clip_imag);

PX_sum = abs(X_sum).^2;

X_DAC_clip = 1/N*fft(s_DAC_clip);
PX_DAC_clip = abs(X_DAC_clip).^2; 

PX_sum(PX_sum==0)=eps;
PX_DAC_clip(PX_DAC_clip==0)=eps;

phases_DAC_clip = angle(X_DAC_clip(ind_ks))'; % 计算DAC clipping后的载波频点相位
phases_DAC_clip = mod(phases_DAC_clip,2*pi); % 相位取到[0, 2*pi]内

err1_phase = abs(phases_DAC_clip - phases)./phases;
disp(err1_phase)
%%
%
figure(1)
for kk=1:Nc
    plot(abs(Xs(:,kk)),'linew',2)
    hold on
end

%
figure(2)
for kk=1:Nc
    subplot(Nc,1,kk)
    plot(real(s_arr(:,kk)),'b-','linew',2)
    hold on
    plot(imag(s_arr(:,kk)),'r-','linew',2)
    xlim([0,N/4])
end

%
figure(3)
plot(real(s_sum),'b-','linew',2)
hold on
plot(imag(s_sum),'r-','linew',2)
xlim([0,N/2])

%
figure(4)
plot(abs(X_sum).^2,'linew',2);

%
figure(5)
subplot(2,1,1)
plot(n,real(s_sum),'b-','linew',2)
hold on
plot(n,imag(s_sum),'r-','linew',2)
xlim([0 N/2])
subplot(2,1,2)
plot(n,real(s_DAC_clip),'b-','linew',2)
hold on
plot(n,imag(s_DAC_clip),'r-','linew',2)
xlim([0 N/2])

%
figure(6)
subplot(2,1,1)
plot(pow2db(PX_sum/P_avg),'linew',1)
subplot(2,1,2)
plot(pow2db(PX_DAC_clip/P_avg),'linew',1)
%%
% 构造脉冲信号
p = zeros(N,1);
p(1) = 1;
X_p = fft(p);
PX_p = X_p.^2;

X_pr = X_p;
X_pr(ind_ks) = 0;
PX_pr = X_pr.^2;
pr = ifft(X_pr);

figure(7)
subplot(2,2,1)
plot(n,p,'m-','linew',2)
subplot(2,2,2)
plot(PX_p,'m-','linew',2)
subplot(2,2,4)
plot(PX_pr,'b-','linew',2)
subplot(2,2,3)
plot(real(pr),'b-','linew',2)
hold on
plot(imag(pr),'r-','linew',2)
%%
s_AG_clip = s_sum;
ft_pr = fft(pr,2*N);
n_iter = 0;
tol = 1e6*eps; % 设定计算误差容限

while(max(max(real(s_AG_clip)),max(imag(s_AG_clip)))-AF >tol || min(min(real(s_AG_clip)),min(imag(s_AG_clip)))+AF< -tol)
    c_real = zeros(N,1);
    c_imag = zeros(N,1);
    ind_real_upper = real(s_AG_clip) > AF;
    ind_real_lower = real(s_AG_clip) < -AF;
    ind_imag_upper = imag(s_AG_clip) > AF;
    ind_imag_lower = imag(s_AG_clip) < -AF;    
    
    c_real(ind_real_upper) = AF - real(s_AG_clip(ind_real_upper));
    c_real(ind_real_lower) = -AF - real(s_AG_clip(ind_real_lower));
    c_imag(ind_imag_upper) = AF - imag(s_AG_clip(ind_imag_upper));
    c_imag(ind_imag_lower) = -AF - imag(s_AG_clip(ind_imag_lower));
    
    c = complex(c_real,c_imag);
    ft_c = fft(c,2*N);
    ft_cc = ft_c.* ft_pr;   

    cc_full = ifft(ft_cc);
    cc = cc_full(1:N);
    
    s_AG_clip = s_AG_clip + cc;

    n_iter = n_iter+1;
    disp(n_iter);
end

X_AG_clip = 1/N*fft(s_AG_clip);
PX_AG_clip = abs(X_AG_clip).^2;

phases_AG_clip = angle(X_AG_clip(ind_ks))'; % 计算AG算法 clipping后的载波频点相位
phases_AG_clip = mod(phases_AG_clip,2*pi); % 相位取到[0, 2*pi]内

err2_phase = abs(phases_AG_clip - phases)./phases;

disp(err1_phase)
disp(err2_phase)
%%
figure(8)
subplot(2,1,1)
plot(real(s_sum),'b-','linew',1)
hold on
plot(imag(s_sum),'r-','linew',1)
subplot(2,1,2)
hold on
plot(real(s_AG_clip),'b-','linew',1)
plot(imag(s_AG_clip),'r-','linew',1)
%
figure(9)
subplot(2,1,1)
plot(pow2db(PX_sum/P_avg),'k-','linew',1)
ylim([-180 0])
subplot(2,1,2)
plot(pow2db(PX_AG_clip/P_avg),'b-','linew',1)
ylim([-180 0])

figure(10)
subplot(2,2,1)
plot(pow2db(PX_sum/P_avg),'k-','linew',1)
subplot(2,2,3)
plot(angle(X_sum),'k-','linew',1)

subplot(2,2,2)
plot(pow2db(PX_AG_clip/P_avg),'b-','linew',1)
% ylim([-180 0])
% xlim([0 200])
subplot(2,2,4)
plot(angle(X_AG_clip),'b-','linew',1)
% ylim([-180 0])
% xlim([0 200])

figure(11)

plot(err1_phase,'m-','linew',2)
hold on
plot(err2_phase,'r-','linew',2)


