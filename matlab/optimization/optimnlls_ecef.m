function F = optimnlls_ecef(X,doa,sensors,combinations)
    
    d = ecef_distance(sensors, X);
    t = zeros(size(combinations,1),1);
    for ii = 1:size(t,1)
       si = combinations(ii,1);
       sj = combinations(ii,2);
       t(ii) = d(si) - d(sj);
    end
    
    F = (doa - t);
end
