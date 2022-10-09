clear;clc;close all;
%%
srate          = 2.4e6; % sampling rate
flms        = 16e3; % the fundemental frequency determined by flux ramp modulation.
n_harmonics = 6; % use n_harmonics to approximate modulated fres.
delay       = 0; % defined by system delay (filter, cryo, etc.) 
mu        = 1e-2; % user defined closed-loop gain
%
lambda = 0.8;
n      = 0:2^14-1;
tn      = n/srate;
% D_squid  = lambda*sin(2*pi*(flms/srate)*n)./(1 + lambda*sin(2*pi*(flms/srate)*n));
D_squid  = lambda*sin(2*pi*flms*tn)./(1 + lambda*sin(2*pi*flms*tn));

%%
% LMS filter
N = 2*n_harmonics+1;
alpha = zeros(1,N);
y = zeros(1,numel(n));

for i = 0:numel(n)-1
    idx = i+1; % array index
%     s_i = [sin(2*pi*flms/srate*i),cos(2*pi*flms/srate*i),...
%         sin(2*pi*flms*2/srate*i),cos(2*pi*flms*2/srate*i),...
%         sin(2*pi*flms*3/srate*i),cos(2*pi*flms*3/srate*i),1];
    s_i = harmonics_gen(n_harmonics,flms,srate,i);
    y(idx) = alpha * s_i';
    e_i = D_squid(idx) - y(idx);
    alpha = alpha + mu*e_i*s_i;   
end
%%
err = D_squid-y;
err_rms = rms(err(2000:end));
err_mean = mean(err(2000:end));
figure(1),clf
subplot(211)
hold on
plot(D_squid,'k-','linew',1)
plot(y,'b.-','linew',1)
xlim([0 1500])
subplot(212)
plot(D_squid-y,'r-','linew',1)
%
