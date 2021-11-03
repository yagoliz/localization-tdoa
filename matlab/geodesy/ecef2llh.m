function llh = ecef2llh(ecef)
% Converts from ECEF to WGS84 lat/lon/height

% Function meat
x = ecef(:,1);
y = ecef(:,2);
z = ecef(:,3);

lon = atan2(y, x);

p = sqrt(x.^2 + y.^2);
th = atan2(geoC.WGS84_A .* z, geoC.WGS84_B .* p);
lat = atan2(z + geoC.wgs84_ep2_b .* sin(th).^3,...
            p - geoC.wgs84_e2_a .* cos(th).^3);

N = geoC.WGS84_A ./ sqrt(1 - geoC.WGS84_ECC_SQ .* sin(lat).^2);
alt = p ./ cos(lat) - N;

llh = [rad2deg(lat), rad2deg(lon), alt];

end

