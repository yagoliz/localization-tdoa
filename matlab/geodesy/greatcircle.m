function d = greatcircle(p0, p1)
    % Returns a great-circle distance in metres between two LLH points,
    % assuming spherical earth_ and _ignoring altitude_. Don't use this if you
    % need a distance accurate to better than 1%

    lat0 = deg2rad(p0(1));
    lon0 = deg2rad(p0(2));
    lat1 = deg2rad(p1(1));
    lon1 = deg2rad(p1(2));
 
    d = geoC.SPHERICAL_R * acos(sin(lat0)*sin(lat1) + ...
        cos(lat0)*cos(lat1)*cos(abs(lon0 - lon1)));   
end