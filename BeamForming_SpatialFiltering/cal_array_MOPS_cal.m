%%
clear;clc;
%
Fs = 250e6; % sampling rate
Nfft = 2048; % FFT-pnts in each PFB.
Ntaps = 4; % taps of one branch in each PFB.
Npfb = 32; % The number of PFB.
Nch = Nfft/2; % channel number.
Npp = 4; % polarization-pair number
Na = 32; % antenna number (each including 2-polarization)
Gc = 4; % complex multiplication factor.
Gr = 1; % real multiplication factor.

ops_rate = Fs/Nfft; % operation rate.(per-sec) 
Mpfb = Npfb*(Ntaps*Nfft*Gr + Nfft/2*log2(Nfft)*Gc); % multipler numbers of PFB
Mcp = Na*Nch*Gc;% multipler numbers of complex phasors after PFB.
Mbf = Nch*Npp*Gc;%  multipler numbers of beamformers after PFBs and phasors.

% mulitpication-operations-per-second under Beamformer Mode.
MOPS_DBF = (Mpfb + Mcp + Mbf)*ops_rate; 

%%
Nx = Na/2*(Na/2-1)/2;
Mcor = Nx*Npp*Nch*Gc;

MOPS_COR = (Mpfb + Mcp + Mcor)*ops_rate; 
