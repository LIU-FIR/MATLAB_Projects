function y = S21_Rp(Qr,Qc,fres,f_offset,phase_dly)
% Compute S21 response when resonator's frequency is offset from ftone 
    S =1-Qr/Qc*1./(1+1i*2*Qr*f_offset./fres); % Ideal S21, fres at real axis.
    y = S.*exp(1i*phase_dly); % uncalibrated S21 measurement

end
