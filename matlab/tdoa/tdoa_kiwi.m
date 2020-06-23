function [doa_samples, iqcorrelate12, corrfactor12, dt] = tdoa_kiwi(iq1,t1,iq2,t2,fs,interpol_factor,corr_type)

if interpol_factor > 1
    iq1 = interp(iq1, interpol_factor);
    iq2 = interp(iq2, interpol_factor);
%     iq1 = resample(iq1, interpol_factor, 1);
%     iq2 = resample(iq2, interpol_factor, 1);
end

j = 10000;
m = 200000;
jump = j * interpol_factor;
maxValue = m * interpol_factor;
iqcorrelate12 = zeros(maxValue/jump,1);
corrfactor12 = zeros(maxValue/jump,1);
dt = zeros(maxValue/jump,1);

data = 1;
for i = 1:jump:maxValue-jump+1
%     [r, lags] = xcorr(iq1(i:i + (jump-1)), iq2(i:i + (jump-1)), round(6371e3*pi/fs), 'coeff');
    [r, lags] = correlate_iq(iq1(i:i + (jump-1)), iq2(i:i + (jump-1)), corr_type, round(6371e3*pi/fs));
    [~,idx] = max(abs(r));
    iqcorrelate12(data) = lags(idx);
    
    timestamps = [t2(i:i+(j-1))'; t1(i:i+(j-1))'];
    corrfactor12(data) = abs(r(idx));
    dt(data) = mean(diff(timestamps));
    data = data + 1;
end
doa_samples = mean(iqcorrelate12);

doa_samples = doa_samples / interpol_factor;
iqcorrelate12 = iqcorrelate12 / interpol_factor;

end

