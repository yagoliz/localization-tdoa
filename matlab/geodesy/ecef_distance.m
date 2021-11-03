function d = ecef_distance(p0, p1)
% Compute cartesian distance from 2 matrices Mx3 where M is the number of
% points

d = sqrt(sum((p0-p1).^2,2));
end