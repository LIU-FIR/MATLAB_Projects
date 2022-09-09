% compute window's function SNR
%

clear; clc;
N = 16384;
hann_win = hann(N);
rect_win = rectwin(N);

enbw_hann = (N*sum(hann_win.^2))./(sum(hann_win)).^2;
enbw_rect = (N*sum(rect_win.^2))./(sum(rect_win)).^2;
disp(['ENBW of hanning window is: ',num2str(enbw_hann)])
disp(['ENBW of rect window is: ',num2str(enbw_rect)])