%%
% This program gives the geometry of IPS-ULA
%

%%
clear; clc;
calarray_xyzs = readmatrix('ljchen_Calarray_MRHC_124_AntennaPosition.csv');
calarray_xs = calarray_xyzs(:,1);
calarray_ys = calarray_xyzs(:,2);
calarray_zs = calarray_xyzs(:,3);
calarray_xys = calarray_xyzs(:,1:2);
%%
figure(1),clf
scatter(calarray_xs,calarray_ys,'bx','linew',1.5)
hold on
plot(0,0,'ko','markerfacecolor','r','markersize',1)
set(gca,'gridlinestyle',':','gridalpha',1)
xlim([-60 60])
ylim([-60 60])
grid on
axis square
title('2-D Calibration Array configuration ')
xlabel('W<—>E (m)')
ylabel('S<—>N (m)')

%%
sarray0_xyzs = calarray_xyzs(1:19,:);
sarray1_xyzs = calarray_xyzs(20:26,:);
sarray2_xyzs = calarray_xyzs(27:33,:);
sarray3_xyzs = calarray_xyzs(34:40,:);
sarray4_xyzs = calarray_xyzs(41:47,:);
sarray5_xyzs = calarray_xyzs(48:54,:);
sarray6_xyzs = calarray_xyzs(55:61,:);
sarray7_xyzs = calarray_xyzs(62:68,:);
sarray8_xyzs = calarray_xyzs(69:75,:);
sarray9_xyzs = calarray_xyzs(76:82,:);
sarray10_xyzs = calarray_xyzs(83:89,:);
sarray11_xyzs = calarray_xyzs(90:96,:);
sarray12_xyzs = calarray_xyzs(97:103,:);
sarray13_xyzs = calarray_xyzs(104:110,:);
sarray14_xyzs = calarray_xyzs(111:117,:);
sarray15_xyzs = calarray_xyzs(118:124,:);

%%
for j = 0:15
    fname = ['MLCA_sarray',num2str(j),'.mat'];
    varname = ['sarray',num2str(j),'_xyzs'];
    save(fname,varname);
end
disp('mat-files are saved')
%%
% valdify via plotting
% save('MLCA_sarray2.mat','sarray2_xyzs');
for j = 0:15
    s1=['sarray',num2str(j),'_xyzs(:,1)'];
    s2=['sarray',num2str(j),'_xyzs(:,2)'];    
    figure(2),clf
    scatter(eval(s1),eval(s2),'bx','linew',1.5)
    title(['sub-array',num2str(j)])
    pause(0.6)
end

%%
% dis = @(x,y) (x.^2 + y.^2).^0.5;
% tol = 1e6*eps;
% % bsxfun(dis,[3,3]',[4,-4]')
% 
% sarray0_xys = zeros(19,2); % xyz-cordinates of sub-array0, the central 19-element array.
% sarray0_xys(1,:) = calarray_xys(1,:); % center of sub-array0
% dis0_xy = calarray_xys - sarray0_xys(1,:);
% dis0 = bsxfun(dis,dis0_xy(:,1),dis0_xy(:,2));
% 
% dis0_th = 11;
% idx_sarr0 = find(dis0 < dis0_th & dis0 > tol); % search other points except for the center.
% 
% sarray0_xys(2:end,:) = calarray_xys(idx_sarr0,:);
% % idx_sarr0_center = find(dis0 < 1e6*eps);
% % sarray0_xys_cen = calarray_xys(idx_sarr0_center,:);
% save('cal_sub_array0.mat','sarray0_xys');
% for kk = 1:19
%     scatter(sarray0_xys(:,1),sarray0_xys(:,2),'ms')
% end
%
