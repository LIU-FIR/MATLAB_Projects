%%
clear;clc;
%%
N = 1024;
n = 0:N-1;
fc = 10e6;
fs = 320e6;
Amax = 1; % clipping level
A = 1.2;

s = A*cos(2*pi*fc/fs*n);
% figure
% plot(n,s,'k-','linew',2)
% % set(gca,'linew',2)
% xlim([0 160])
%
s_clip = s;
s_clip(s_clip>Amax)=Amax;
s_clip(s_clip<-Amax)=-Amax;
figure(1)
plot(n,s,'k-','linew',2)
hold on
plot(n,s_clip,'b-','linew',2)
xlim([0 160])
%
clip_waveform = s_clip-s;
figure(2)
plot(n,clip_waveform,'r-','linew',2)
xlim([0 160])

%
sx = 1/N*fft(s);
psx = abs(sx).^2;
sx_clip = 1/N*fft(s_clip);
psx_clip = abs(sx_clip).^2; 
sx_clip_waveform = 1/N*fft(clip_waveform);
psx_clip_waveform = abs(sx_clip_waveform).^2;

figure(3)
subplot(3,1,1)
plot(pow2db(psx))
subplot(3,1,2)
plot(pow2db(psx_clip))
subplot(3,1,3)
plot(pow2db(psx_clip_waveform))
%
pulse_t = zeros(1,N);
pulse_t(1) = 1;
sx_pulse_t = fft(pulse_t);
psx_pulse_t = sx_pulse_t.^2;
ind_k = 33;

sx_pulse_t_rev = sx_pulse_t;
sx_pulse_t_rev(ind_k) = 0;
psx_pulse_t_rev = sx_pulse_t_rev.^2;
pulse_t_rev = real(ifft(sx_pulse_t_rev));

figure(4)
subplot(2,2,1)
plot(n,pulse_t,'m-','linew',2)
subplot(2,2,2)
plot(psx_pulse_t,'m-','linew',2)
subplot(2,2,4)
plot(psx_pulse_t_rev,'m-','linew',2)
subplot(2,2,3)
plot(pulse_t_rev,'m-','linew',2)
%%
s_tmp = s;
c = zeros(1,length(s));
ft_pulse_t_rev = fft(pulse_t_rev,2*N);
n_iter = 0;

% for kk = 1:1
% c = zeros(1,length(s));
% ind_upper = s_tmp > Amax;
% ind_lower = s_tmp < -Amax;
% c(ind_upper) = Amax - s_tmp(ind_upper);
% c(ind_lower) = -Amax - s_tmp(ind_lower);
% 
% ft_c = fft(c,2*N);
% ft_cc = ft_c.* ft_pulse_t_rev;   
% 
% cc_full = real(ifft(ft_cc));
% cc = cc_full(1:N);
% % cc = conv(c,pulse_t_rev,'same');
% s_tmp = s_tmp + cc;
% % c = zeros(1,length(s));
% disp(max(cc))
% end
% %%
% figure(7)
% plot(cc,'r-','linew',2)
while(max(s_tmp)> Amax || min(s_tmp)< -Amax)
    c = zeros(1,length(s));
    ind_upper = s_tmp > Amax;
    ind_lower = s_tmp < -Amax;
    c(ind_upper) = Amax - s_tmp(ind_upper);
    c(ind_lower) = -Amax - s_tmp(ind_lower);

    ft_c = fft(c,2*N);
    ft_cc = ft_c.* ft_pulse_t_rev;   

    cc_full = real(ifft(ft_cc));
    cc = cc_full(1:N);
    % cc = conv(c,pulse_t_rev,'same');
    s_tmp = s_tmp + cc;
    % c = zeros(1,length(s));
    disp(max(cc))
end
%%
% s1 = s + cc;
figure(5)
plot(s,'k-','linew',2)
hold on
plot(s_tmp,'b-','linew',2)
%%
sx_tmp = 1/N*fft(s_tmp);
psx_tmp = abs(sx_tmp).^2;

figure(6)
subplot(2,1,1)
plot(pow2db(psx),'k-','linew',2)
subplot(2,1,2)
plot(pow2db(psx_tmp),'b-','linew',2)

%%
% ind = [33 97 161 225 289 353 417 481 545 609 673 737 801 865 929 993];
% S = sx;
% Sp = sx_clip-sx_clip_waveform;
% 
% PS = sum(abs(S).^2);
% PSp = sum(abs(Sp).^2);
% 
% disp(PS)
% disp(PSp)
