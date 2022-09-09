
clear;clc;

%%
s1 = 2000;
s2 = 4000;
xs = -100:.1:100;
ys = -100:.1:100;
A_x = exp(-xs.^2/s1); 
A_y = exp(-ys.^2/s2); 
A = A_x'*A_y;
%%

figure(1)
colormap('jet')
imagesc(A)
%%
figure(2),clf
contour(A,7)
xlim([1000 2000])
ylim([1000 2000])
xticks([])
yticks([])
xlabel('Pipeline granularity')
ylabel('Parallel computation')