function y = harmonics_gen(Nord,f,srate,ti)
    f_Hs = f*(1:Nord);
    y_sin = upsample(sin(2*pi*f_Hs/srate*ti),2);
    y_cos = upsample(cos(2*pi*f_Hs/srate*ti),2,1);
    y = [y_sin+y_cos 1];
end
