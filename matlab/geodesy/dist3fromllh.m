function d = dist3fromllh(llh0, llh1)
% Compute cartesian distance from (latitude, longitude, altitude)
p0 = llh2ecef(llh0);
p1 = llh2ecef(llh1);

d = ecef_distance(p0, p1);
end