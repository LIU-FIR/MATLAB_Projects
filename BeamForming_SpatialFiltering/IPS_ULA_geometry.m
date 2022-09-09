%%
% This program gives the geometry of IPS-ULA
%

%%
clear; clc;
N_Elements = 592; % number of array elements in each cylinder telescope.
N_Cylinders = 3; % number of cylindrical telescopes.

fA = 327e6; % radio frequency, (Hz)
fB = 654e6;
cspeed = 3e8; % light speed, in m/s
lambdaA = cspeed/fA; % radio wavelength, (unit: m)
lambdaB = cspeed/fB;
%
d_NS = lambdaB/2; % The spacing of ULA elements.
d_EW = 40;
IPS_ULA_Geometry_NS = (-N_Elements/2:1:N_Elements/2-1)'.*[d_NS d_NS d_NS];
IPS_ULA_Geometry_EW = ones(N_Elements,1).*[-d_EW 0 d_EW];
% m = 8; % The number of array elements.
% K = 74; % 
IPS_ULA_xs = IPS_ULA_Geometry_EW(:);
IPS_ULA_ys = IPS_ULA_Geometry_NS(:);
%%
figure(1)
scatter(IPS_ULA_xs,IPS_ULA_ys,'bx')
hold on
plot(0,0,'ko','markerfacecolor','r','markersize',5)
set(gca,'gridlinestyle',':','gridalpha',1)
xlim([-50 50])
grid on
axis square
title(['IPS ULA configuration: N-S spacing \sim ',num2str(round(d_NS,2)),'m'])
xlabel('W<—>E (m)')
ylabel('S<—>N (m)')

%%

a_rs = exp((0:m-1)'.*(-1i*2*pi*d_NS*sin(deg2rad(angle_rs))/lambdaA)); % direction vector along radio source's direction
array_vec = exp((0:m-1)'*(-1i*2*pi*d_NS*sin(deg2rad(angles))/lambdaA)); % array direction-vector 
X = a_rs'*array_vec; % 8-element synthesized pattern response over all directions.
X_max = max(abs(X));
X_db = mag2db(abs(X)./X_max);% synthesized power pattern resnponse 


%%
figure(1),clf
subplot(211),hold on
plot(angles,X_db,'k-','linew',2)
ylim([-50 2])
yrange = get(gca,'ylim');
p1 = plot([angles(tdx1), angles(tdx1)],[yrange(1),yrange(2)],'r--','linew',1);
p2 = plot([angles(tdx2), angles(tdx2)],[yrange(1),yrange(2)],'r--','linew',1);
legend([p1],{'-3dB points'})
title('8-element synthesized power pattern @327MHz')
xlabel('angle ($\theta$): degree','Interpreter','latex')
ylabel('power pattern: dB')
grid on
xticks([angle_lim(1),angle_null_1,angle_3dB_down_1,0,angle_3dB_down_2,angle_null_2,angle_lim(2)])
subplot(212),hold on
plot(angles,X_db,'k-','linew',2)
plot(angles,X1_db,'b-.','linew',1)
plot(angles,X2_db,'m-.','linew',1)
grid on
ylim([-50 2])
title('3dB-width shifted power patterns @327MHz')
xlabel('angle ($\theta$): degree','Interpreter','latex')
ylabel('power pattern: dB')
