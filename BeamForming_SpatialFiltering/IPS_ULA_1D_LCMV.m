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
%
% LCMV digital beamformer: In -3dB synthesized power pattern range, set one
% direction of source and one interfere direction.
%
K = 74; % The number of digital array elements
L = 200; % signal time-length (taps).
srate = 250e6; % sampling rate.
snr = 10; % signal-to-noise ratio
inr = 10; % interference-to-noise ratio

Nc = 401; % DBF computed directions.
angle_3dB_pnt = [rad2deg(angles_zl(tdx1)) rad2deg(angles_zl(tdx2))];
angles_dbf_span = linspace(angle_3dB_pnt(1),angle_3dB_pnt(2),Nc);

d_eff = m*d; % effective delay between each synthesized pattern pair.
% angles_rs = angles_dbf_span(floor(Nc/2)+1); % source direction
angles_rs = rad2deg(angle_zl_rs);% source direction
idx1_inf = floor(Nc/2)+1-100;
idx2_inf = floor(Nc/2)+1+100;
angles_inf = [angles_dbf_span(idx1_inf),angles_dbf_span(idx2_inf)]; % interference directions
%%
v_rs = exp(d_eff.*(0:K-1)'*1i*2*pi*sin(deg2rad(angles_rs)));
v_inf = exp(d_eff.*(0:K-1)'*1i*2*pi*sin(deg2rad(angles_inf)));
t = (0:L-1)./srate; % signal time interval.

x_rs = sqrt(10^(snr/10))*v_rs*exp(1i*2*pi*fA*t); %complex signal has sqrt(2) gain.
x_inf = sqrt(10^(snr/10)/2)*v_inf*(randn(length(angles_inf),L)+1i*randn(length(angles_inf),L));
noise = (randn(K,L)+1i*randn(K,L))/sqrt(2); % noise of elements.

Z = x_inf + noise;

R = Z*Z'; % covariance matrix

w_lcmv = (R\v_rs)/(v_rs'*(R\v_rs)); % LCMV-DBF weights
w_fix = v_rs; % Fixed-DBF weights

v_span = exp(d_eff.*(0:K-1)'*1i*2*pi*sin(deg2rad(angles_dbf_span)));

B = abs(w_lcmv'*v_span);
B_db = mag2db(B./max(B,[],2));

B1 = abs(w_fix'*v_span);
B1_db = mag2db(B1./max(B1,[],2));
%%
figure(3),clf
hold on
% plot(angles_dbf,y1_db,'r','linew',1)
% plot(angles_dbf,y2_db,'b','linew',1)
% plot(angles_dbf,Y_rs_db(9,:));
xlim([10 50])
ylim([-70 2])
p1=plot(angles_dbf_span,B1_db,'b-','linew',1);
p2=plot(angles_dbf_span,B_db,'m-','linew',1);
plot(rad2deg(angles_zl),Y_zl_db,'k-','linew',1.5)
% plot(angles,SP_db,'k-','linew',1.5)
p3 = plot(angles_dbf_span(idx1_inf),-70,'rx','linew',2);
plot(angles_dbf_span(idx2_inf),-70,'rx','linew',2)
yrange = get(gca,'ylim');
plot([rad2deg(angles_zl(tdx1)), rad2deg(angles_zl(tdx1))],[yrange(1),yrange(2)],'k--','linew',1)
plot([rad2deg(angles_zl(tdx2)), rad2deg(angles_zl(tdx2))],[yrange(1),yrange(2)],'k--','linew',1)
xticks([angles_zl_lim(1),angle_3dB_down_1,round(angles_dbf_span(idx1_inf),2),...
    rad2deg(angle_zl_rs),round(angles_dbf_span(idx2_inf),2),angle_3dB_down_2,angles_zl_lim(2)])
% xticklabels()
grid on
legend([p1 p2 p3],{'Fixed-DBF','LCMV-DBF','RFI locations'})

xlabel('angle ($\theta$): degree','Interpreter','latex')
ylabel('power pattern: dB')
title('IPS digital beamformer: LCMV-DBF v.s. fixed-DBF')
%%
