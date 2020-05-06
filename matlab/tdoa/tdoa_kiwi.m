function [doa_samples] = tdoa_kiwi(signal1,signal2,interpol_factor,corr_type)

if interpol_factor > 1
    signal1corr = interp(signal1, interpol_factor);
    signal2corr = interp(signal2, interpol_factor);
   
else
    signal1corr = signal1;
    signal2corr = signal2;
end

% Compute the lags
[corr_signal, lags] = correlate_iq(signal1, signal2, corr_type, 1);
[corr_signalinterp, lagsinterp] = correlate_iq(signal1corr, signal2corr, corr_type, 1);

% DOA (without interpolation)
[corr_val, idx] = max(corr_signal);
doa = lags(idx);

% DOA (with interpolation)
[corr_val_interp, idx_interp] = max(corr_signalinterp);
doa_samples = lagsinterp(idx_interp)/interpol_factor;

end

