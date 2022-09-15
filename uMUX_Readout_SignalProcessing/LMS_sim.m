clear;clc;close all;

Fs          = 2.4e6;
Flms        = 16e3;
n_harmonics = 3;
delay       = 22; % defined by system delay (filter, cryo, etc.) 
gain        = 1e-2;
%%
% Integrator 
ol          = tf([1], [1 -1], 1/Fs);

% Add harmonic stages in parallel
for i = 1:n_harmonics
    w0    = 2*pi*i*Flms/Fs;
    phi0  = angle(exp(1j*w0));
    ol    = ol + tf([cos(w0*(1 + delay)), -cos(phi0*delay)], ...
	        [1, -2*cos(w0), 1, zeros(1, delay)], 1/Fs);
end

% Closed loop transfer function
cl = feedback(gain*ol, 1);

% Simulate SQUID
lambda = 0.8;
n      = 0:2^14-1;
t      = n/Fs;
squid  = lambda*sin(2*pi*(Flms/Fs)*n)./(1 + lambda*sin(2*pi*(Flms/Fs)*n));

%
% figure(1)
% bode(cl)
%
y = lsim(cl, squid, t);
%%
figure(2),clf
plot(squid,'k-','linew',1)
hold on
plot(y,'b-','linew',1)
xlim([0 1200])