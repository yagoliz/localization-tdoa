function [ iq_corr, lags ] = correlate_iq( iq1, iq2, corr_strategy, max_lag)
%correlate_iq correlates two complex iq signals according to the specified strategy
% iq1: first complex IQ signal
% iq2: second complex IQ signal
% corr_strategy ={'abs', 'dphase'}
%   'abs': use absolute value of iq signals
%   'dphase': use differential phase of iq signals
% smoothing_factor

    switch corr_strategy
        
        case 'abs'
            abs1 = remove_mean(abs(iq1));
            abs2 = remove_mean(abs(iq2));

            [abs_corr, lags] = xcorr(abs1, abs2, 'coeff');
            
%             abs_corr = abs_corr ./ max(abs(abs_corr)); %normalize
            iq_corr = abs_corr;
                        
        case 'dphase'
            d_phase1 = diff(unwrap(angle(iq1)));
            d_phase2 = diff(unwrap(angle(iq2)));
            
            d_phase1 = [ 0; d_phase1(1:length(d_phase1)) ];
            d_phase2 = [ 0; d_phase2(1:length(d_phase2)) ];

            % We need to flatten the d_phase shape
            d_phase1 = detrend(d_phase1);
            d_phase2 = detrend(d_phase2);
            
            [d_phase_corr, lags] = xcorr(d_phase1, d_phase2, 'coeff');
            
%             d_phase_corr = d_phase_corr ./ max(abs(d_phase_corr)); %normalize
            iq_corr = d_phase_corr;
            
        case 'phase'
            d_phase1 = detrend(unwrap(angle(iq1)));
            d_phase2 = detrend(unwrap(angle(iq2)));
           

            d_phase1 = remove_mean(d_phase1);
            d_phase2 = remove_mean(d_phase2);        
            
            [d_phase_corr, lags] = xcorr(d_phase1, d_phase2, 'coeff');
            
%             d_phase_corr = d_phase_corr ./ max(abs(d_phase_corr)); %normalize
            iq_corr = d_phase_corr;
            
        case 'iq'
            [iq_corr, lags] = xcorr(remove_mean(iq1), remove_mean(iq2), 'coeff');
            
%             iq_corr = iq_corr ./ max(abs(iq_corr));
            
        otherwise
            error('Correlation strategy not supported (choose abs, dphase, phase or iq)');
    end 

end

