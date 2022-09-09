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
% angle_span = 45; % max sky angle span,(degree)
theta_rs = 5; % radio source(rs) direction, with respect to zenith.(degree)
phi_rs = 1; % radio source(rs) direction, with respect to x-axis.(degree)
angle_res = rad2deg(lambdaA/(m-1)/d); % array aperture resolution

% angles = -45:1:45; % scan the angle
angle_lim = [-45 45];
angles = linspace(angle_lim(1),angle_lim(2),K);
% angles_sky = -angle_span:angle_res:angle_span;
angle_rs = angles(K/2+1);
angle_delta = diff(angles(1:2));
angle_rs1 = angle_rs + angle_delta;
angle_rs2 = angle_rs - angle_delta;

array_vec = exp((0:m-1)'*(1i*2*pi*d*sin(deg2rad(angles))/lambdaA)); % direction vector of radio source
%%
w_rs = exp(1i*2*pi*d*sin(deg2rad(angle_rs))/lambdaA.*(0:m-1)'); % weighted vector along radio source's direction
w_rs1 = exp(1i*2*pi*d*sin(deg2rad(angle_rs1))/lambdaA.*(0:m-1)'); % weighted vector along radio source's direction
w_rs2 = exp(1i*2*pi*d*sin(deg2rad(angle_rs2))/lambdaA.*(0:m-1)'); % weighted vector along radio source's direction

y = w_rs'*array_vec; % response
y1 = w_rs1'*array_vec; % response
y2 = w_rs2'*array_vec; % response

% synthesized response matrix: each col is associated with one direction.
bmp_mat = array_vec'*array_vec; 
figure(1),clf
plot(angles,abs(bmp_mat),'linew',1)
title('Synthesized beam patterns')
grid on

bmp_syn_db = mag2db(abs(y/max(abs(y))));
bmp1_syn_db = mag2db(abs(y1/max(abs(y1))));
bmp2_syn_db = mag2db(abs(y2/max(abs(y2))));
tdx1 = dsearchn(bmp_syn_db(1:K/2)',-3); % beam pattern -3dB point
tdx2 = dsearchn(bmp_syn_db(K/2+1:end)',-3)+K/2;% beam pattern -3dB point
tdx3 = dsearchn(bmp_syn_db(1:K/2)',-0.5); % beam pattern -0.5dB point
tdx4 = dsearchn(bmp_syn_db(K/2+1:end)',-0.5)+K/2;% beam pattern -0.5dB point
xlabel('angle ($\theta$): degree','Interpreter','latex')
ylabel('power pattern: dB')

figure(2),clf
subplot(211),hold on
plot(angles,bmp_syn_db,'k-','linew',2)
plot(angles,bmp1_syn_db,'b--','linew',2)
plot(angles,bmp2_syn_db,'m--','linew',2)
grid on
xlabel('angle ($\theta$): degree','Interpreter','latex')
ylabel('power pattern: dB')
title('A single  synthesized beam pattern')
subplot(212),hold on
plot(angles,bmp_syn_db,'k-','linew',2)
yrange = get(gca,'ylim');
p1 = plot([angles(tdx1), angles(tdx1)],[yrange(1),yrange(2)],'r--','linew',1);
p2 = plot([angles(tdx2), angles(tdx2)],[yrange(1),yrange(2)],'r--','linew',1);
p3 = plot([angles(tdx3), angles(tdx3)],[yrange(1),yrange(2)],'k--','linew',1);
p4 = plot([angles(tdx4), angles(tdx4)],[yrange(1),yrange(2)],'k--','linew',1);
title('A single synthesized beam pattern')
xticks([angle_lim(1),round(angles(tdx1),2),round(angles(tdx3),2),0,...
    round(angles(tdx4),2),round(angles(tdx2),2),angle_lim(2)])
grid on
disp(['-3dB power pattern width: ',num2str(abs(angles(tdx1)-angles(tdx2)))])
disp(['-0.5dB power pattern width: ',num2str(abs(angles(tdx3)-angles(tdx4)))])
xlabel('angle ($\theta$): degree','Interpreter','latex')
ylabel('power pattern: dB')
legend([p1 p3],{'-3dB points','-0.5dB points'})
%%
% digital beamformer: 16 beams in -1dB angle range.
J = 256; % computed directions.
angle_1dB_pnt = [angles(tdx3) angles(tdx4)];
angles_dbf = linspace(angle_1dB_pnt(1),angle_1dB_pnt(2),J);
d_syn = m*d;

W_rs = exp((0:K-1)'*(1i*2*pi*d_syn*sin(deg2rad(angles_dbf(2:16:end))))/lambdaA);

array_syn_vec = exp((0:K-1)'*(1i*2*pi*d_syn*sin(deg2rad(angles_dbf))/lambdaA)); % direction vector of synthesized array
w_syn_rs1 = exp(1i*2*pi*d_syn*sin(deg2rad(angles_dbf(J/2+1)))/lambdaA.*(0:K-1)');
w_syn_rs2 = exp(1i*2*pi*d_syn*sin(deg2rad(angles_dbf(J/2+1+16)))/lambdaA.*(0:K-1)');
% w_syn_rs = ones(74,1);
Z_rs = W_rs'*array_syn_vec;
Z_rs_db = mag2db(abs(Z_rs)./max(abs(Z_rs),[],2));

z1 = w_syn_rs1'*array_syn_vec;
z2 = w_syn_rs2'*array_syn_vec;
y1_db = mag2db(abs(z1)./max(abs(z1)));
y2_db = mag2db(abs(z2)./max(abs(z2)));
%%
figure(3),clf
hold on
% plot(angles_dbf,y1_db,'r','linew',1)
% plot(angles_dbf,y2_db,'b','linew',1)
% plot(angles_dbf,Y_rs_db(9,:));
plot(angles_dbf,Z_rs_db','linew',1);
plot(angles,bmp_syn_db,'k-','linew',2)
yrange = get(gca,'ylim');
% plot([angles(tdx1), angles(tdx1)],[yrange(1),yrange(2)],'r--','linew',1)
% plot([angles(tdx2), angles(tdx2)],[yrange(1),yrange(2)],'r--','linew',1)
plot([angles(tdx3), angles(tdx3)],[yrange(1),yrange(2)],'k--','linew',1)
plot([angles(tdx4), angles(tdx4)],[yrange(1),yrange(2)],'k--','linew',1)
% xticks([angle_lim(1),round(angles(tdx1),2),round(angles(tdx3),2),...
%     0,round(angles(tdx4),2),round(angles(tdx2),2),angle_lim(2)])
xticks([angle_lim(1),round(angles(tdx3),2),...
    0,round(angles(tdx4),2),angle_lim(2)])
% xticklabels()
grid on
xlim([-6 6])
ylim([-70 2])
xlabel('angle ($\theta$): degree','Interpreter','latex')
ylabel('power pattern: dB')
title('16 digital beams,(rect-windowed weights)')
%%
hann_win = hann(K);
wh_syn_rs1 = w_syn_rs1.*hann_win;
W_hann_rs = W_rs.*hann_win;
z1_hann = wh_syn_rs1'*array_syn_vec;
z1_hann_db = mag2db(abs(z1_hann)./max(abs(z1_hann)));

Y_hann_rs = W_hann_rs'*array_syn_vec;
Y_hann_rs_db = mag2db(abs(Y_hann_rs)./max(abs(Y_hann_rs),[],2));

figure(4),clf
hold on

plot(angles_dbf,Y_hann_rs_db','linew',1);
plot(angles,bmp_syn_db,'k-','linew',2)
yrange = get(gca,'ylim');

plot([angles(tdx3), angles(tdx3)],[yrange(1),yrange(2)],'k--','linew',1)
plot([angles(tdx4), angles(tdx4)],[yrange(1),yrange(2)],'k--','linew',1)
% xticks([angle_lim(1),round(angles(tdx1),2),round(angles(tdx3),2),...
%     0,round(angles(tdx4),2),round(angles(tdx2),2),angle_lim(2)])
xticks([angle_lim(1),round(angles(tdx3),2),...
    0,round(angles(tdx4),2),angle_lim(2)])
% xticklabels()
grid on
xlim([-6 6])
ylim([-70 2])
xlabel('angle ($\theta$): degree','Interpreter','latex')
ylabel('power pattern: dB')
title('16 digital beams,(hann-windowed weights)')
%%
% using DFT to compute pattern.
dsita = deg2rad(diff(angles_dbf(1:2)));
sita_res = lambdaA/(K*d_syn);
Nx = ceil(sita_res/dsita*K);

x = exp(1i*2*pi*d_syn*sin(deg2rad(angles_dbf(J/2)))/lambdaA.*(0:K-1)');
sita_ap = rad2deg(lambdaA/(K*d_syn))*((0:K-1)-K/2);
X = fft(x,Nx);
Y = abs(X)./max(abs(X));
Y_db = mag2db(Y);
sita = dsita*(0:Nx-1);
figure(5),clf
plot(rad2deg(sita),fftshift(Y_db),'k-','linew',1)
