clear;clc;

filename = 'c:\udp-pipeline\real_data\sat_samples_mini.bin';

fd = fopen(filename);
data = fread(fd,'int8');
fclose(fd);
disp('Reading done');
%%

clear;clc;
filename = 'c:\udp-pipeline\real_data\power_spec8192FT_4097ch_128000fra.mat';
mf = matfile(filename);

tmp = mf.X_power_frames(:,600:2000);
%%
figure(1),clf
imagesc(tmp)
axis xy
colormap('jet')

%%
N = 8192;

data = reshape(data,N,length(data)/N);

[~, nseg] = size(data);
%%
X_power_frames = abs(fft(data)).^2;
disp('fft is done')


%%
X_power_frames = X_power_frames(1:4097,:);
save('c:\udp-pipeline\real_data\power_spec8192FT_4097ch_128000fra.mat','X_power_frames','-v7.3');
disp('saving done')
%%
figure(1),clf
imagesc(X_power_frames)
axis xy
colormap('jet')
% xlim([600 1200])
%%
X_power_avg = 0;
for k = 1:nseg
    X = fft(data(:,k));
    X_power_avg = abs(X).^2 + X_power_avg;
end
%
disp('computing done')
%%
X_power_avg_dB = 20*log10(X_power_avg);


%%
fax = linspace(1024,1536,N/2+1);
figure(1),clf

%%
load('c:\Users\FIR.LIU\udp-xiaoyun_notebook\data_fft8192.mat')
load('c:\udp-pipeline\real_data\data_fft8192_linear.mat')
%%
error = X_power_avg_dB(1:N/2+1)'- Pavg_8192;
error1 = X_power_avg(1:N/2+1)'- Pavg_8192_linear;
%%
figure(2),clf
subplot(211)
plot(X_power_avg_dB(1:N/2+1),'b-','linew',1)
subplot(212)
plot(fax, error1./Pavg_8192_linear,'k-','linew',1)