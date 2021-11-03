function F = optimnlls_llh(X,doa,sensors,combinations,keep)
    lat = X(1);
    lon = X(2);
    
    d = zeros(size(sensors,1),1);
    for ii = 1:size(sensors,1)
        d(ii) = deg2km(distance(lat, lon, sensors(ii,1), sensors(ii,2)));
    end
    
    t = zeros(size(combinations,1),1);
    for ii = 1:size(t,1)
        si = combinations(ii,1);
        sj = combinations(ii,2);
        t(ii) = (d(si) - d(sj))*1000;
    end
    
    F = (doa(keep) - t(keep));
end