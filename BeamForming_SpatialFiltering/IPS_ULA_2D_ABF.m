%%
% IPS ULA 8-to-1 analog synthesized beam pattern (2-D)
% by Liu Fei
% 2022,April

%%
clear; clc;

fA = 327e6; % radio frequency, (Hz)
fB = 654e6;
cspeed = 3e8; % light speed, in m/s
lambdaA = cspeed/fA; % radio wavelength, (unit: m)
lambdaB = cspeed/fB;

%
m = 8; % The number of array elements.
d = lambdaB/2; % The spacing of elements.
Nd_ze = 101; % 
Nd_az = 101;

% set (zenith_angle,azimuth_angle) of  the radio source
arad_ze_rs = deg2rad(30); % zenith_angle: zl
arad_az_rs = deg2rad(180);% azimuth_angle: az

adeg_ze_lim = [0 60];
adeg_az_lim = [150 210];

adegs_ze = linspace(adeg_ze_lim(1),adeg_ze_lim(2),Nd_ze);
adegs_az = linspace(adeg_az_lim(1),adeg_az_lim(2),Nd_az);
arads_ze = deg2rad(linspace(adeg_ze_lim(1),adeg_ze_lim(2),Nd_ze));
arads_az = deg2rad(linspace(adeg_az_lim(1),adeg_az_lim(2),Nd_az));
% angles_sky = -angle_span:angle_res:angle_span;

% a_x = zeros(Nd_az,Nd_ze); 
a_x = zeros(Nd_ze,Nd_az); 

for i_az = 1:Nd_az
    for i_ze = 1:Nd_ze
        a_x(i_ze,i_az)=cos(arads_az(i_az))*sin(arads_ze(i_ze));
    end
end
av_x = a_x(:)'; % x-component direction vector

A_x = exp(d*(0:m-1)'*(1i*2*pi*av_x./lambdaA)); % x-component array direction matrix 
w_rs_x = exp(d*(0:m-1)'*1i*2*pi*sin(arad_ze_rs)*cos(arad_az_rs)/lambdaA); % weighted vector along radio source's direction

Y_x = abs(w_rs_x'*A_x); % response
Y_x_2d = reshape(Y_x,Nd_ze,Nd_az);

%%

%%
% Use imagesc(), add some -3dB region.

Y_x_2d_max = max(Y_x_2d,[],'all');
Y_3db_down = Y_x_2d_max*10^(-3/20); 

figure(1),clf
hold on
colormap('jet'),colorbar
imagesc(Y_x_2d), % turn i-axis to azimuth-angle, and j-axis to zenith-angle
axis xy % turn i-j to x-y, Up and right increasing.
title('IPS ULA 8-element synthesized beam pattern @327MHz')
xlabel('azimuth angle: degree')
ylabel('zenith angle: degree')
% xticks(linspace(angles_az_lim(1),angles_az_lim(2),Nd_az))
xticks([1 51 101])
xticklabels({'150','180','210'})
yticks([1 51 101])
yticklabels({'0','30','60'})
% plot([51 51],[10 20],'k--','linew',1)
xlim([1 101])
ylim([1 101])
cRange = caxis;
contour(Y_x_2d,[Y_3db_down Y_3db_down],'k--')
caxis(cRange);
% clim([min(Y_x_2d_dB,[],'all') max(Y_x_2d_dB,[],'all')])
%%
% [phi_grid,theta_grid] = meshgrid(linspace(angle_al_lim(1),angle_al_lim(2),Nd));
[ze_grid,az_grid] = meshgrid(linspace(adeg_ze_lim(1),adeg_ze_lim(2),Nd_ze),...
                                linspace(adeg_az_lim(1),adeg_az_lim(2),Nd_az));
figure(2),clf
% subplot(131)
hold on
colormap('jet'),colorbar
surf(ze_grid,az_grid,Y_x_2d)
title('IPS ULA 8-element synthesized beam pattern @327MHz')
xlabel('azimuth angle: degree')
ylabel('zenith angle: degree')
