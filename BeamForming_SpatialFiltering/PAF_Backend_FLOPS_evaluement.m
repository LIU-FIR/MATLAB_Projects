% Evaluate FLOPS in PAF Backend system
% FLOPS = Number of multipliers * operatation_rate
% 
% 2022/05/20 Liu Fei

%%

clear;clc;

% dfeq_ZoomFFT = 0.1e3;
fs = 500e6;
os_ratio_str = '4/3';
os_ratio = str2num(os_ratio_str);


M = 2048; % 1-st PFB-channels
R = 4;
K = 512; % Zoom FFT-length,



f_clk1 = fs/M*os_ratio;
f_clk2 = fs/M/K*os_ratio;

dt_CoarseFFT = 1/f_clk1;

dfeq_ZoomFFT = fs/M/K*os_ratio;



N_ele = 300; % array physical elements.
N_pol_PerElement = 2; % polarizations for each element.
N_inputs = N_ele * N_pol_PerElement; % all inputs from array elements.
N_RealTime_Beams = 90; % Number of real-time digital beams (single polarization.)
N_CoarseChans = M/2 ;
N_FineChans_PerCoarseChan = K; % complex DBF signal input.

N_Mul_CoarsePFB_PerInput = M/2*log2(M) + M*R; % Multipliers of one single coarse PFB.
N_Mul_DBF_PerCoarseChan_PerBeam = 2*N_ele; % Multipliers of one single beam.(2 polarizations)
N_Mul_CoarseSptr_PerCoarseChan_PerBeam = 4;
N_Mul_Covar_PerCoarseChan = 4*N_ele*N_ele; % Multipliers of one covariance block(4 matrices).
N_Mul_ZoomFFT_PerCoarseChan_PerBeam = K/2*log2(K);
N_Mul_ZoomSptr_PerFineChan_PerCoarseChan_PerBeam = 4;

%

N_Mul_CoarsePFBs = N_Mul_CoarsePFB_PerInput * N_inputs *4;

N_Mul_DBFs = N_Mul_DBF_PerCoarseChan_PerBeam... 
            * N_CoarseChans...
            * N_RealTime_Beams *4;

N_Mul_CoarseSptrs = N_Mul_CoarseSptr_PerCoarseChan_PerBeam... 
                    * N_CoarseChans...
                    * N_RealTime_Beams *4;

N_Mul_ZoomFFT = 2 * N_Mul_ZoomFFT_PerCoarseChan_PerBeam...
                  * N_CoarseChans * N_RealTime_Beams *4;
                
N_Mul_ZoomSptrs = N_Mul_ZoomSptr_PerFineChan_PerCoarseChan_PerBeam...
                * N_FineChans_PerCoarseChan...
                * N_CoarseChans...
                * N_RealTime_Beams * 4;

N_Mul_Covar = N_Mul_Covar_PerCoarseChan * N_CoarseChans *4;
            
%%
flops = f_clk1 * N_Mul_CoarsePFBs + f_clk1*N_Mul_DBFs + f_clk1*N_Mul_Covar...
    + f_clk1 * N_Mul_CoarseSptrs + f_clk2 * N_Mul_ZoomFFT + f_clk2 * N_Mul_ZoomSptrs;

flops_str = sprintf('FLOPs is: %2.4e ',flops);

clc;
disp(['M = ',num2str(M),';K = ',num2str(K),...
    ';Fine spectral resolution(kHz):',num2str(dfeq_ZoomFFT/1e3)])
disp(['M = ',num2str(M),';os_ratio = ',os_ratio_str,...
    ';CoarseFFT time resolution(us):',num2str(dt_CoarseFFT*1e6)])
disp(flops_str)