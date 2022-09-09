%%
% 计算Fs=10 kHz,Nfft=2048, 501点长度的Hanning window的ENBW
clear;clc;
L = 501;
Nfft = 2048;
fs = 10e3;
fres = fs/Nfft;

fcn_hanning = hann(L);
fcn_rect = rectwin(L);

enbw_bin_hann = enbw(hann(L))*Nfft/L; % 计算Nfft点的 ENBW和L点的ENBW之间存在比例关系

enbw_hz_hann = enbw_bin_hann*fs/Nfft;

X_mag_hanning = abs(fft(fcn_hanning,Nfft));
X_mag_rect = abs(fft(fcn_rect,Nfft));
X_logmag_norm_hanning = 20*log10(X_mag_hanning/max(X_mag_hanning));
X_logmag_norm_rect = 20*log10(X_mag_rect/max(X_mag_rect));
%
fk = -fs/2:fres:fs/2-fres;
figure(1)
h = plot(fk,fftshift(X_logmag_norm_hanning),'k-');
hold on
h1 = plot(fk,fftshift(X_logmag_norm_rect),'b-');
x_p = [-enbw_hz_hann/2 -enbw_hz_hann/2 enbw_hz_hann/2 enbw_hz_hann/2];
y_limit = get(gca,'ylim');
y_p = [y_limit(1), y_limit(2), y_limit(2), y_limit(1)];
patch(x_p,y_p,'k','facealpha',.25,'EdgeColor','#D95319','linew',2)
set(gca,'xlim',[-100 100],'ylim',[-80 5])
set(h,'linew',2)
set(h1,'linew',2)
grid on
xlabel('Frequency: Hz','FontSize',14','Interpreter','latex')
ylabel('$\mathrm{Magnitude: 20log_{10}|W(k)/|W(0)|}$','Interpreter','latex','FontSize',14)
title({'Equivalent noise bandwidth'},{ '$\mathrm{Hanning Window length\, L = 501, N_{FFT} = 2048, F_s = 10\, kHz}$'},'Interpreter','latex','FontSize',14)
text(x_p(3),-58,'$\bf \leftarrow ENBW = 30 Hz$','FontSize',14,'Interpreter','latex')
text(x_p(1),-58,'$\bf \rightarrow$','FontSize',14,'HorizontalAlignment','right','Interpreter','latex')
% hold on
%
%%

%%