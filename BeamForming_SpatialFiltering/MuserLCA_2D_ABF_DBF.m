%%
% MUSER-LCA
% MUSER-Low-frequency-calibration-array analog synthesized beam pattern (2-D)
% 19-element sub-array(station)
% by Liu Fei
% 2022,May

%%
% 19-element array in the center
% indexing: sub-array0
clear; clc;
load('MLCA_sarray0.mat')
load('MLCA_sarray1.mat')
load('MLCA_sarray2.mat')
load('MLCA_sarray3.mat')
load('MLCA_sarray4.mat')
load('MLCA_sarray5.mat')
load('MLCA_sarray6.mat')
load('MLCA_sarray7.mat')
load('MLCA_sarray8.mat')
load('MLCA_sarray9.mat')
load('MLCA_sarray10.mat')
load('MLCA_sarray11.mat')
load('MLCA_sarray12.mat')
load('MLCA_sarray13.mat')
load('MLCA_sarray14.mat')
load('MLCA_sarray15.mat')
%
% freqx is 50MHz,100MHz,150MHz,200MHz,250MHz,300MHz,350MHz,400MHz
fAs = 1e6*(50:50:400);
cspeed = 3e8; % light speed, in m/s
idx_array = input('Choose subarray[0-15]:');
% idx_array = 0;
%%
for tt = 2:2
fA = fAs(tt); % radio frequency, (Hz)
lambdaA = cspeed/fA; % radio wavelength, (unit: m)

sc1=['sarray',num2str(idx_array),'_xyzs(1,1)'];
sc2=['sarray',num2str(idx_array),'_xyzs(1,2)']; 
s1=['sarray',num2str(idx_array),'_xyzs(:,1)'];
s2=['sarray',num2str(idx_array),'_xyzs(:,2)'];  
xs = eval(s1)-eval(sc1);
ys = eval(s2)-eval(sc2);
% xs = sarray0_xyzs(:,1)-sarray0_xyzs(1,1);
% ys = sarray0_xyzs(:,2)-sarray0_xyzs(1,2);
% xs = sarray1_xyzs(:,1)-sarray1_xyzs(1,1);
% ys = sarray1_xyzs(:,2)-sarray1_xyzs(1,2);
% xs = sarray7_xyzs(:,1)-sarray7_xyzs(1,1);
% ys = sarray7_xyzs(:,2)-sarray7_xyzs(1,2);
%
%
Nd_ze = 201; % sampling points in zenith-angle direction to compute.
Nd_az = 201; % sampling points in azimuth-angle direction to compute.

arad_ze_rs = deg2rad(30); % radio source(rs) direction relative to zenith.(rad)
arad_az_rs = deg2rad(180); % radio source(rs) direction relative to x-axis.(rad)

adeg_ze_lim = [0 60]; % altitude angle range, (degree)
adeg_az_lim = [150 210]; % azimuth angle range, (degree)

% angles sampling points in range (rad)
adegs_ze=linspace(adeg_ze_lim(1),adeg_ze_lim(2),Nd_ze);
adegs_az=linspace(adeg_az_lim(1),adeg_az_lim(2),Nd_az);
adeg_res_ze = diff(adegs_ze(1:2));
adeg_res_az = diff(adegs_az(1:2));
arads_ze = deg2rad(linspace(adeg_ze_lim(1),adeg_ze_lim(2),Nd_ze)); 
arads_az = deg2rad(linspace(adeg_az_lim(1),adeg_az_lim(2),Nd_az));

% direction factor
% a_x = zeros(Nd_az,Nd_ze);
a_x = zeros(Nd_ze,Nd_az);
% a_y = zeros(Nd_az,Nd_ze);
a_y = zeros(Nd_ze,Nd_az);

% compute the direction factors of all (Nd_az * Nd_al) points 
for i_az = 1:Nd_az
    for i_ze = 1:Nd_ze
%         a_x(i_az,i_ze)=cos(arads_az(i_az))*sin(arads_ze(i_ze));
%         a_y(i_az,i_ze)=sin(arads_az(i_az))*sin(arads_ze(i_ze));
        a_x(i_ze,i_az)=cos(arads_az(i_az))*sin(arads_ze(i_ze));
        a_y(i_ze,i_az)=sin(arads_az(i_az))*sin(arads_ze(i_ze));        
    end
end
aV_x = a_x(:)'; % x-component of direction factors in vector.
aV_y = a_y(:)'; % y-component of direction factors in vector.

A_x = exp(xs*(1i*2*pi*aV_x./lambdaA)); % x-component matrix of all directions
A_y = exp(ys*(1i*2*pi*aV_y./lambdaA)); % y-component matrix of all directions
A = A_x.*A_y; % The composed matrix of all directions.
%
% the weight vector associated with the radio source direction sensored by 2D planar array.
w_rs_x = exp(xs*1i*2*pi*cos(arad_az_rs)*sin(arad_ze_rs)/lambdaA);
w_rs_y = exp(ys*1i*2*pi*sin(arad_az_rs)*sin(arad_ze_rs)/lambdaA);
w_rs = w_rs_x.*w_rs_y;
%
Y = w_rs'*A; % response in all directions
Y_x = w_rs'*A_x;
Y_y = w_rs'*A_y;

% Y_2d = abs(reshape(Y,Nd_az,Nd_ze));
% Y_x_2d = abs(reshape(Y_x,Nd_ze,Nd_az));
% Y_y_2d = abs(reshape(Y_y,Nd_ze,Nd_az));
Y_2d = abs(reshape(Y,Nd_ze,Nd_az));
%

% Y_2d = Y_2d'; % turn i-axis to zenith-angle; j-axis to azimuth-angle.
Y_2d_max = max(Y_2d,[],'all'); % compute -3dB value.
Y_3db_down = Y_2d_max*10^(-3/20); 
%
figure(tt),clf
hold on
colormap('jet'),colorbar
title(['Response of sub-array-',num2str(idx_array),' of MUSER-LCA @',num2str(fA/1e6),'MHz'])
% axis square
xlabel('azimuth angle: degree')
ylabel('zenith angle: degree')
imagesc(Y_2d)
axis xy
xticks([1 floor(Nd_az/2)+1 Nd_az])
xticklabels({'150','180','210'})
yticks([1 floor(Nd_az/2)+1 Nd_ze])
yticklabels({'0','30','60'})
xlim([1 Nd_az])
ylim([1 Nd_ze])
cRange = caxis;
M0=contour(Y_2d,[Y_3db_down Y_3db_down],'k--','linew',1);
caxis(cRange);

end
%%
% Compute -3dB angle-range
idx_centr = [floor(Nd_ze/2)+1,floor(Nd_az/2)+1]; % central-point of image's index.
M0q = round(M0); % fit -3dB cord to int-index.
M0q_x = M0q(1,:);M0q_y = M0q(2,:);
idx_3db_x = find(M0q_x==idx_centr(1));
idx_3db_y = find(M0q_y==idx_centr(2));
range_3db_x = abs(M0q(1,idx_3db_y(1))-M0q(1,idx_3db_y(2)));
range_3db_y = abs(M0q(2,idx_3db_x(1))-M0q(2,idx_3db_x(2)));
disp('3dB angle-range (degree)')
disp([range_3db_x*adeg_res_az,range_3db_y*adeg_res_ze])
% [ze_grid,az_grid] = meshgrid(linspace(angles_ze_lim(1),angles_ze_lim(2),Nd_ze),...
%                                 linspace(angles_az_lim(1),angles_az_lim(2),Nd_az));
% subplot(122),hold on
% colormap('jet'),colorbar
% surf(az_grid,ze_grid,Y_2d)
% title(['Response of sub-array-',num2str(idx_array),' of MUSER-LCA @100MHz'])
% axis xy
% xlabel('azimuth angle: degree')
% ylabel('zenith angle: degree')
% v = [18,18];
% contourf(az_grid,zl_grid,mag_Y_2d,v,'k--','linew',2)
%%
figure(3),clf
subplot(121),hold on
scatter(eval(s1),eval(s2),'bx','linew',1.5)
grid on
% set(gca,'xlim',[-15 15],'ylim',[-15 15])
axis square
xlabel('W<—>E (m)')
ylabel('S<—>N (m)')
title(['MUSER-LCA planar sub-array',num2str(idx_array),' configuration'])

subplot(122),hold on
colormap('jet'),colorbar
title(['Response of sub-array-',num2str(idx_array),' of MUSER-LCA @',num2str(fA/1e6),'MHz'])
% axis square
xlabel('azimuth angle: degree')
ylabel('zenith angle: degree')
imagesc(Y_2d)
axis xy
xticks([1 floor(Nd_az/2)+1 Nd_az])
xticklabels({'150','180','210'})
yticks([1 floor(Nd_az/2)+1 Nd_ze])
yticklabels({'0','30','60'})
xlim([1 Nd_az])
ylim([1 Nd_ze])
cRange = caxis;
contour(Y_2d,[Y_3db_down Y_3db_down],'k--','linew',1)
caxis(cRange);

%% 
% 2) compute 16 digital beams.

xys_subarr_centr = zeros(16,2);
% all 16 subarray phase centers coordinates
for ks = 0:15
    
    scx=['sarray',num2str(ks),'_xyzs(1,1)'];
    scy=['sarray',num2str(ks),'_xyzs(1,2)']; 
    xys_subarr_centr(ks+1,:) = [eval(scx),eval(scy)];
end

xs_sc = xys_subarr_centr(:,1);
ys_sc = xys_subarr_centr(:,2);

Ndbf_ze = 201; % sampling points in zenith-angle direction to compute.
Ndbf_az = 201; % sampling points in azimuth-angle direction to compute.

arad_ze_rs = deg2rad(30); % radio source(rs) direction relative to zenith.(rad)
arad_az_rs = deg2rad(180); % radio source(rs) direction relative to x-axis.(rad)

adeg_ze_bf_lim = [0 60]; % altitude angle range, (degree)
adeg_az_bf_lim = [150 210]; % azimuth angle range, (degree)

% angles sampling points in range (rad)
adegs_ze_bf = linspace(adeg_ze_bf_lim(1),adeg_ze_bf_lim(2),Ndbf_ze);
adegs_az_bf = linspace(adeg_az_bf_lim(1),adeg_az_bf_lim(2),Ndbf_az);

arads_ze_bf = deg2rad(adegs_ze_bf); 
arads_az_bf = deg2rad(adegs_az_bf);

% direction factor
a_bf_x = zeros(Ndbf_az,Ndbf_ze);
a_bf_y = zeros(Ndbf_az,Ndbf_ze);

% compute the direction factors of all (Nd_az * Nd_al) points 
for i_az = 1:Ndbf_az
    for i_ze = 1:Ndbf_ze
        a_bf_x(i_az,i_ze)=cos(arads_az_bf(i_az))*sin(arads_ze_bf(i_ze));
        a_bf_y(i_az,i_ze)=sin(arads_az_bf(i_az))*sin(arads_ze_bf(i_ze));
    end
end
aV_bf_x = a_bf_x(:)'; % x-component of direction factors in vector.
aV_bf_y = a_bf_y(:)'; % y-component of direction factors in vector.

A_bf_x = exp(xs_sc*(1i*2*pi*aV_bf_x./lambdaA)); % x-component matrix of all directions
A_bf_y = exp(ys_sc*(1i*2*pi*aV_bf_y./lambdaA)); % y-component matrix of all directions
A_bf = A_bf_x.*A_bf_y; % The composed matrix of all directions.
%
% the weight vector associated with the radio source direction sensored by 2D planar array.
w_bf_rs_x = exp(xs_sc*1i*2*pi*cos(arad_az_rs)*sin(arad_ze_rs)/lambdaA);
w_bf_rs_y = exp(ys_sc*1i*2*pi*sin(arad_az_rs)*sin(arad_ze_rs)/lambdaA);
w_bf_rs = w_bf_rs_x.*w_bf_rs_y;
%
Y_bf = w_bf_rs'*A_bf; % response in all directions
Y_bf_x = w_bf_rs'*A_bf_x;
Y_bf_y = w_bf_rs'*A_bf_y;

Y_bf_2d = abs(reshape(Y_bf,Ndbf_az,Ndbf_ze));
Y_bf_x_2d = abs(reshape(Y_bf_x,Ndbf_az,Ndbf_ze));
Y_bf_y_2d = abs(reshape(Y_bf_y,Ndbf_az,Ndbf_ze));
%
Y_bf_2d = Y_bf_2d'; % turn i-axis to zenith-angle; j-axis to azimuth-angle.
Y_bf_2d_max = max(Y_bf_2d,[],'all');
Y_bf_3db_down = Y_bf_2d_max*10^(-3/20);
%%
figure(4)
hold on
colormap('jet'),colorbar
title('Response of digital beamformer of MUSER-LCA @100MHz')
% axis square
xlabel('azimuth angle: degree')
ylabel('zenith angle: degree')
imagesc(Y_bf_2d)
axis xy
xticks([1 floor(Ndbf_az/2)+1 Ndbf_az])
% xticklabels({'165','180','195'})
xticklabels({'150','180','210'})
yticks([1 floor(Ndbf_ze/2)+1 Ndbf_ze])
% yticklabels({'165','180','195'})
yticklabels({'0','30','60'})

xlim([1 Ndbf_az])
ylim([1 Ndbf_ze])
cRange = caxis;
% contour(Y_bf_2d',[Y_bf_3db_down Y_bf_3db_down],'k--','linew',1)
M1 = contour(Y_2d,[Y_3db_down Y_3db_down],'k--','linew',1);
caxis(cRange);
M2 = contour(Y_bf_2d,[Y_bf_3db_down Y_bf_3db_down],'k--','linew',1.5);
caxis(cRange);
%% 
% compute and plot multiple digital beams 
% 
% 1) compute 3dB range indices
arad_res_ze_bf = diff(arads_ze_bf(1:2)); % angle resolution of sampling points.
arad_res_az_bf = diff(arads_az_bf(1:2));

idx_centr = [floor(Ndbf_ze/2)+1,floor(Ndbf_az/2)+1]; % central-point of image's index.

M1q = round(M1); % fit -3dB cord to int-index.
M1q_x = M1q(1,:);M1q_y = M1q(2,:);
idx_3db_x = find(M1q_x==idx_centr(1));
idx_3db_y = find(M1q_y==idx_centr(2));
range_3db_x = abs(M1q(1,idx_3db_y(1))-M1q(1,idx_3db_y(2)));
range_3db_y = abs(M1q(2,idx_3db_x(1))-M1q(2,idx_3db_x(2)));


M2q = round(M2);
M2q_x = M2q(1,:);M2q_y = M2q(2,:);
idx_dbf_3db_x = find(M2q_x==idx_centr(1));
idx_dbf_3db_y = find(M2q_y==idx_centr(2));
range_dbf_3db_x = abs(M2q(1,idx_dbf_3db_y(1))-M2q(1,idx_dbf_3db_y(2)));
range_dbf_3db_y = abs(M2q(2,idx_dbf_3db_x(1))-M2q(2,idx_dbf_3db_x(2)));
%
% Non-overlapped digtal beams.
% N_dbeam_x = ceil(range_3db_x/range_dbf_3db_x);
N_dbeam_y = ceil(range_3db_y/range_dbf_3db_y);
N_dbeam_x = 3;
N_dbeam = N_dbeam_x * N_dbeam_y;


%
Ys_bf_2d = zeros(Ndbf_az,Ndbf_ze,N_dbeam);

[arad_ze_grid,arad_az_grid] = meshgrid(arad_ze_rs + (-floor(N_dbeam_x/2):floor(N_dbeam_x/2))*range_dbf_3db_x*arad_res_ze_bf,...
                            arad_az_rs + (-floor(N_dbeam_y/2):floor(N_dbeam_y/2))*range_dbf_3db_y*arad_res_az_bf);
arad_list = [arad_ze_grid(:),arad_az_grid(:)];

%%
% compute direction-factor-matrix
ak_bf_x = zeros(Ndbf_az,Ndbf_ze);
ak_bf_y = zeros(Ndbf_az,Ndbf_ze);

% compute the direction factors of all (Nd_az * Nd_al) points 
for i_az = 1:Ndbf_az
    for i_ze = 1:Ndbf_ze
        ak_bf_x(i_az,i_ze)=cos(arads_az_bf(i_az))*sin(arads_ze_bf(i_ze));
        ak_bf_y(i_az,i_ze)=sin(arads_az_bf(i_az))*sin(arads_ze_bf(i_ze));
    end
end
aVk_bf_x = ak_bf_x(:)'; % x-component of direction factors in vector.
aVk_bf_y = ak_bf_y(:)'; % y-component of direction factors in vector.

Ak_bf_x = exp(xs_sc*(1i*2*pi*aVk_bf_x./lambdaA)); % x-component matrix of all directions
Ak_bf_y = exp(ys_sc*(1i*2*pi*aVk_bf_y./lambdaA)); % y-component matrix of all directions
Ak_bf = Ak_bf_x.*Ak_bf_y; % The composed matrix of all directions.

%%
for nn = 1:N_dbeam
    
    arad_ze_now = arad_list(nn,1);
    arad_az_now = arad_list(nn,2);
    %
    % the weight vector associated with the radio source direction sensored by 2D planar array.
    wk_bf_rs_x = exp(xs_sc*1i*2*pi*cos(arad_az_now)*sin(arad_ze_now)/lambdaA);
    wk_bf_rs_y = exp(ys_sc*1i*2*pi*sin(arad_az_now)*sin(arad_ze_now)/lambdaA);
    wk_bf_rs = wk_bf_rs_x.*wk_bf_rs_y;
    %
    Yk_bf = wk_bf_rs'*Ak_bf; % response in all directions
    Yk_bf_2d = abs(reshape(Yk_bf,Ndbf_az,Ndbf_ze));
    %
    Yk_bf_2d = Yk_bf_2d'; % turn i-axis to zenith-angle; j-axis to azimuth-angle.
    Ys_bf_2d(:,:,nn) = Yk_bf_2d;

end
disp('All digital beams computation over.')

span3db_x = idx_centr(2)-range_3db_x/2:idx_centr(2)+range_3db_x/2;
span3db_y = idx_centr(1)-range_3db_y/2:idx_centr(1)+range_3db_y/2;

%%
figure(5),clf
hold on
colormap('jet'),colorbar
title('Response of digital beamformer of MUSER-LCA @100MHz')
% axis square
xlabel('azimuth angle: degree')
ylabel('zenith angle: degree')

% imagesc(sum(Ys_bf_2d(:,:,5),3))
% imagesc(Ys_bf_2d(:,:,2))
xticks([1 floor(Ndbf_az/2)+1 Ndbf_az])
% xticklabels({'165','180','195'})
xticklabels({'150','180','210'})
yticks([1 floor(Ndbf_ze/2)+1 Ndbf_ze])
% yticklabels({'165','180','195'})
yticklabels({'0','30','60'})

% write first frame of .gif
imagesc(Ys_bf_2d(:,:,1))
cRange = caxis;
contour(Y_2d,[Y_3db_down Y_3db_down],'k--','linew',1);
caxis(cRange);
gif('test.gif')
% write the rest frames of .gif
for nn = 2:N_dbeam
    imagesc(Ys_bf_2d(:,:,nn))
    cRange = caxis;
    contour(Y_2d,[Y_3db_down Y_3db_down],'k--','linew',1);
    caxis(cRange);
    gif('DelayTime',0.2,'LoopCount',1)
end
%

%%
% Ym_bf_2d = ones(Ndbf_az,Ndbf_ze);
% 
% [idx_xg,idx_yg] = meshgrid((-1:1)*range_dbf_3db_x, (-1:1)*range_dbf_3db_y);
% idx_grid = [idx_xg(:)+idx_centr(2),idx_yg(:)+idx_centr(1)];
% % idx_centr(1) % ze
% % idx_centr(2) % az
% for jj = 1:9
%     spanj_x = idx_grid(jj,1)-range_dbf_3db_x/2:idx_grid(jj,1)+range_dbf_3db_x/2;
%     spanj_y = idx_grid(jj,2)-range_dbf_3db_x/2:idx_grid(jj,2)+range_dbf_3db_x/2;
%     Ym_bf_2d(spanj_x,spanj_y) = Ys_bf_2d(spanj_x,spanj_y,jj) + Ym_bf_2d(spanj_x,spanj_y);
% end
% %%
% span3db_x = idx_centr(2)-range_3db_x/2:idx_centr(2)+range_3db_x/2;
% span3db_y = idx_centr(1)-range_3db_y/2:idx_centr(1)+range_3db_y/2;
% figure(6),clf
% % caxis(cRange);
% ax1 = axes;
% im1 = imagesc(Ys_bf_2d(span3db_x,span3db_y,9),cRange);
% im1.AlphaData = 0.5;
% axis xy
% figure(7),clf
% ax2 = axes;
% im2 = imagesc(Ys_bf_2d(span3db_x,span3db_y,8),cRange);
% im2.AlphaData = 0.5;
% axis xy
% linkaxes([ax1,ax2])
% colormap('jet'),colorbar
