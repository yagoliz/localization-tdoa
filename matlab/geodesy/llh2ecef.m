function coord = llh2ecef(llh)
% Converts from WGS84 lat/lon/height to ellipsoid-earth ECEF

lat = deg2rad(llh(:,1));
lng = deg2rad(llh(:,2));
alt = llh(:,3);

slat = sin(lat);
slng = sin(lng);
clat = cos(lat);
clng = cos(lng);

d = sqrt(1 - (slat .* slat .* geoC.WGS84_ECC_SQ));
rn = geoC.WGS84_A ./ d;

x = (rn + alt) .* clat .* clng;
y = (rn + alt) .* clat .* slng;
z = (rn .* (1 - geoC.WGS84_ECC_SQ) + alt) .* slat;

coord = [x, y, z];
end