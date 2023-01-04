function y = Eta_gen(S21_ms,fres,f,rotation_mode)
%   S21_ms: S21 measurement,normally fres is not located on the imagniary
%   axis.
%   rotation mode: 'real',when fres is rotated to real axis;'imag',when fres is rotated to imaginary axis
    fzc = f-fres; % move original fres to zero frequency,
    idx_fres = find(abs(S21_ms)==min(abs(S21_ms))); % find "shifted" resonance frequency.
    Eta = (fzc(idx_fres+1)-fzc(idx_fres-1))./(S21_ms(idx_fres+1)-S21_ms(idx_fres-1)); % compute Eta by RF (SMuRF)
    if isequal(rotation_mode,'imag')
        y = Eta;
    elseif isequal(rotation_mode,'real')
        y = Eta.*1i;
    else
        disp('illegal rotation_mode');
        y = nan;
    end
end
