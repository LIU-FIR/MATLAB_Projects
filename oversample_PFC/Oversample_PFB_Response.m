%%
% Polyphase Filter Bank: PFB的Matlab实现
%%
clear; clc;

M = 8;
R = 20;
P = 40;

L = M*R*P;
N = M*R-1;
Fc = 1/M; % 截止频率=pi/32
% Fp = 1/16;
% Fst = 1/32;
% Ap = 1;
% Astop = 40;
Hf = fdesign.lowpass('N,Fc',N,Fc);
% Hf = fdesign.lowpass('Fp,Fst,Ap,Ast');
Hd2 = design(Hf,'window','window',{@chebwin,50}, ...
            'systemobject',true);
% hfvt = fvtool(Hd2,'Color','White');
h = Hd2.Numerator;
%%
%画出最大抽取和过采样的滤波器响应
figure(1)
plot((-0.5:1/2048:.5-1/2048)*M,fftshift(20*log10(abs(fft(h,2048)))),'r')
grid on

%%

B = 2*pi/M;
w0 = 2*B;
% w1 = pi/64+3*B;
w1 = w0+B/64;
w2 = w0-B/64;
% x = exp(1i*w0*(0:L-1))+1.5*exp(1i*w1*(0:L-1))+4*exp(1i*w2*(0:L-1));
x = exp(1i*w0*(0:L-1))+exp(-1i*w0*(0:L-1));
% x=ones(1,L);

y_fir = pfb_fir(x,h,M,R);
y_os_fir = os2_pfb_fir(x,h,M,R);
disp('end here.')
%%
ini_num1 = size(y_fir,2);
ini_num2 = size(y_os_fir,2);
y_psd_mean = pfb_psd_mean(y_fir,ini_num1);
y_os_psd_mean = pfb_psd_mean(y_os_fir,ini_num2);
% y_os_psd_mean1 = mean(abs(conj(fft(y_os_fir))).^2,2);
figure(2)
stem(y_psd_mean,'r*')
hold on
stem(y_os_psd_mean,'b')

%%
% PFB的频率响应：在全带宽[0,2*pi]生成多个复正弦信号，以扫频方式画出每个正弦信号的频率响应。
fspan = linspace(0,M*B,8*101);
t=(0:M*R*P-1);
chan = zeros(M,length(fspan));
chan_os = zeros(M,length(fspan));
% chan_d = zeros
for k = 1:length(fspan)
    z = exp(1i*fspan(k)*t);
%     z2 = exp(1i*fspan(k)/2*t);
    z_fir = pfb_fir(z,h,M,R);
    z2_os_fir = os2_pfb_fir(z,h,M,R);
    init1 = size(z_fir,2);
    init2 = size(z2_os_fir,2);
    z_psd = pfb_psd_mean(z_fir,init1);
    z2_os_psd = pfb_psd_mean_os(z2_os_fir,init2);
    for p = 1:M
        chan(p,k) = z_psd(p);
        chan_os(p,k) = z2_os_psd(p);
    end
end
disp('end response computation')
%%
figure(3)
for kk = 1:2
    plot(fspan/pi, db(chan(kk,:)))
    hold on
end
% plot(fspan/pi, db(chan(2,:)))
xlim([0,fspan(end)/pi])
ylim([-150 1])
xlabel('Normalized frequency: (\times \pi /rad/sample)')
ylabel('Channel Power: dB')
grid on
%%
figure(4)
for kk = 1:5
    plot(fspan/pi, db(chan_os(kk,:)))
    hold on    
    pause
end
% plot(fspan/pi/2, db(chan_os(2,:)))
hold on
% plot(fspan/pi, db(chan(2,:)),'r-')
xlim([0,fspan(end)/pi])
ylim([-150 1])
xlabel('Normalized frequency: (\times \pi /rad/sample)')
ylabel('Channel Power: dB')
grid on

%%
% fspan = linspace(0,B,101);
% t=(0:M*R*P-1);
% chan_rp = zeros(1,length(fspan));
% chan_os_rp = zeros(1,length(fspan));
% z = zeros(1,M*R*P);
% for k = 1:length(fspan)
%     z = exp(1i*fspan(k)*t);
%     z_fir = pfb_fir(z,h,M,R);
%     z_os_fir = oversampled_pfb_fir(z,h,M,R);
%     init1 = size(z_fir,2);
%     init2 = size(z_os_fir,2);
%     z_psd = pfb_psd_mean(z_fir,init1); 
%     z_os_psd = pfb_psd_mean(z_os_fir,init2); 
%     chan_rp(k)=z_psd(1);
%     chan_os_rp(k)= z_os_psd(1);
% end
% disp('end response');
% %%
% figure(3)
% plot(fspan/pi, db(chan_rp))
% hold on
% plot(fspan/pi, db(chan_os_rp),'r-')
% ylim([-150 1])
% xlabel('Normalized frequency: (\times \pi /rad/sample)')
% ylabel('Channel Power: dB')
% grid on