%%
% This program gives a brief model description of narrow-band beamforming.
% Part 1: Source & Sky
% Part 2: Sampling (spatial and temporal)
% Part 3: Receiving and spatial filtering

%%
clear; clc;

fA = 327e6; % radio frequency, (Hz)
fB = 654e6;
cspeed = 3e8; % light speed, in m/s
lambdaA = cspeed/fA; % radio wavelength, (unit: m)
lambdaB = cspeed/fB;
%%
m = 8; % The number of array elements.
d = lambdaB/2; % The spacing of elements.
K = 74; % 
Nd_zl = 401; % The directions of array.
Nd_az = 11; % The directions of array.
% angle_span = 45; % max sky angle span,(degree)
angle_zl_rs = deg2rad(30); % zenith_angle: zl
angle_az_rs = deg2rad(180);% azimuth_angle: az
% angle_res = rad2deg(lambdaA/(m-1)/d); % array aperture resolution

% angles = -45:1:45; % scan the angle
angles_zl_lim = [0 60];
angles_zl = deg2rad(linspace(angles_zl_lim(1),angles_zl_lim(2),Nd_zl));
angles_az_lim = [150 210];
angles_az = deg2rad(linspace(angles_az_lim(1),angles_az_lim(2),Nd_az));
% angles_sky = -angle_span:angle_res:angle_span;
% angle_rs = angles(floor(Nd_zl/2)+1); % from zenith.
% angle_delta = diff(angles(1:2));
%%
a_x = zeros(Nd_az,Nd_zl); 
% a_y = zeros(Nd_az,Nd_al);
for i_az = 1:Nd_az
    for i_zl = 1:Nd_zl
        a_x(i_az,i_zl)=cos(angles_az(i_az))*sin(angles_zl(i_zl));
%         a_y(i_az,i_al)=sin(angles_az(i_az))*sin(angles_al(i_al));
    end
end
av_x = a_x(:)'; % x-component direction vector

A_x = exp(d*(0:m-1)'*(1i*2*pi*av_x./lambdaA)); % x-component array direction matrix 

w_rs_x = exp(d*(0:m-1)'*1i*2*pi*sin(angle_zl_rs)*cos(angle_az_rs)/lambdaA); % weighted vector along radio source's direction

Y_x = abs(w_rs_x'*A_x); % response

Y_x_2d = reshape(Y_x,Nd_az,Nd_zl);

idx_az_rs = floor(Nd_az/2)+1;

Y_zl = abs(Y_x_2d(idx_az_rs,:));
Y_zl_max = max(abs(Y_zl));
Y_zl_db = mag2db(Y_zl./Y_zl_max);

%%
tdx1 = dsearchn(Y_zl_db(1:floor(Nd_zl/2))',-3); % power pattern -3dB point
tdx2 = dsearchn(Y_zl_db(floor(Nd_zl/2)+1:end)',-3)+floor(Nd_zl/2);% power pattern -3dB point

angle_3dB_down_1 = round(rad2deg(angles_zl(tdx1)),2);
angle_3dB_down_2 = round(rad2deg(angles_zl(tdx2)),2);

figure(1),clf
hold on
plot(rad2deg(angles_zl),Y_zl_db,'k-','linew',2)
ylim([-50 2])
yrange = get(gca,'ylim');
p1 = plot([rad2deg(angles_zl(tdx1)), rad2deg(angles_zl(tdx1))],[yrange(1),yrange(2)],'r--','linew',1);
p2 = plot([rad2deg(angles_zl(tdx2)), rad2deg(angles_zl(tdx2))],[yrange(1),yrange(2)],'r--','linew',1);
legend([p1],{'-3dB points'})
xticks([angles_zl_lim(1),angle_3dB_down_1,rad2deg(angle_zl_rs),angle_3dB_down_2,angles_zl_lim(2)])
title('IPS 8-element synthesized power pattern (zenith-angle scan) @327MHz')
xlabel('Zenith angle: (degree)')
ylabel('Power pattern: (dB)')
grid on


%%
% synthesized response matrix: each col is associated with one direction.
% bmp_mat = array_vec'*array_vec; 
%
% digital beamformer: 16 beams in -3dB angle range.
Nc = 401; % DBF computed directions.
angle_3dB_pnt = [rad2deg(angles_zl(tdx1)) rad2deg(angles_zl(tdx2))];
angles_dbf_span = linspace(angle_3dB_pnt(1),angle_3dB_pnt(2),Nc);
d_eff = m*d; % effective delay between each synthesized pattern pair.

% presumed direction of 16 radio sources
angles_rs_16 = angles_dbf_span(1:25:end); 
angles_rs_16 = angles_rs_16(1:16);
%%
Wd = exp((0:K-1)'*(1i*2*pi*d_eff*sin(deg2rad(angles_rs_16)))/lambdaA);

array_dbf_vec = exp((0:K-1)'*(1i*2*pi*d_eff*sin(deg2rad(angles_dbf_span))/lambdaA)); % direction vector of synthesized array
% w_syn_rs1 = exp(1i*2*pi*d_eff*sin(deg2rad(angles_dbf(Nc/2+1)))/lambdaA.*(0:K-1)');
% w_syn_rs2 = exp(1i*2*pi*d_eff*sin(deg2rad(angles_dbf(Nc/2+1+16)))/lambdaA.*(0:K-1)');
% w_syn_rs = ones(74,1);
Y_rs = Wd'*array_dbf_vec;
Y_rs_db = mag2db(abs(Y_rs)./max(abs(Y_rs),[],2));
% 
%%
figure(3),clf
hold on
% plot(angles_dbf,y1_db,'r','linew',1)
% plot(angles_dbf,y2_db,'b','linew',1)
% plot(angles_dbf,Y_rs_db(9,:));
plot(angles_dbf_span,Y_rs_db','linew',1);
plot(rad2deg(angles_zl),Y_zl_db,'k-','linew',2)
yrange = get(gca,'ylim');
plot([rad2deg(angles_zl(tdx1)), rad2deg(angles_zl(tdx1))],[yrange(1),yrange(2)],'r--','linew',1)
plot([rad2deg(angles_zl(tdx2)), rad2deg(angles_zl(tdx2))],[yrange(1),yrange(2)],'r--','linew',1)
xticks([angles_zl_lim(1),angle_3dB_down_1,rad2deg(angle_zl_rs),angle_3dB_down_2,angles_zl_lim(2)])
% xticks([angle_lim(1),round(angles_zl(tdx1),2),...
%     0,round(angles(tdx2),2),angle_lim(2)])
grid on
xlim([10 50])
ylim([-70 2])
xlabel('angle ($\theta$): degree','Interpreter','latex')
ylabel('power pattern: dB')
title('16 digital beams,(rect-windowed weights)')
%%
hann_win = hann(K);
Wd_hann = Wd.*hann_win;


Y_hann_rs = Wd_hann'*array_dbf_vec;
% Y_hann_rs_db = mag2db(abs(Y_hann_rs)./max(abs(Y_hann_rs),[],2));
Y_hann_rs_db = mag2db(abs(Y_hann_rs)./max(abs(Y_rs),[],2));

figure(4),clf
hold on
plot(angles_dbf_span,Y_hann_rs_db','linew',1);
plot(rad2deg(angles_zl),Y_zl_db,'k-','linew',2)
yrange = get(gca,'ylim');

plot([rad2deg(angles_zl(tdx1)), rad2deg(angles_zl(tdx1))],[yrange(1),yrange(2)],'r--','linew',1)
plot([rad2deg(angles_zl(tdx2)), rad2deg(angles_zl(tdx2))],[yrange(1),yrange(2)],'r--','linew',1)
% xticks([angle_lim(1),round(angles(tdx1),2),round(angles(tdx3),2),...
%     0,round(angles(tdx4),2),round(angles(tdx2),2),angle_lim(2)])
xticks([angles_zl_lim(1),angle_3dB_down_1,rad2deg(angle_zl_rs),angle_3dB_down_2,angles_zl_lim(2)])
grid on
xlim([10 50])
ylim([-70 2])
xlabel('angle ($\theta$): degree','Interpreter','latex')
ylabel('power pattern: dB')
title('16 digital beams,(hann-windowed weights)')
%%
% using DFT to compute pattern.
% dsita = deg2rad(diff(angles_dbf_span(1:2)));
% sita_res = lambdaA/(K*d_eff);
% Nx = ceil(sita_res/dsita*K);
% 
% X = exp(1i*2*pi*d_eff*sin(deg2rad(angles_dbf_span(Nc/2)))/lambdaA.*(0:K-1)');
% sita_ap = rad2deg(lambdaA/(K*d_eff))*((0:K-1)-K/2);
% X = fft(X,Nx);
% Y = abs(X)./max(abs(X));
% Y_db = mag2db(Y);
% sita = dsita*(0:Nx-1);
% figure(5),clf
% plot(rad2deg(sita),fftshift(Y_db),'k-','linew',1)
