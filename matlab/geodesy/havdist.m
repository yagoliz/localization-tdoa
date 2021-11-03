function d = havdist(ll0, ll1)
% Compute haversine distance (up to 1% precision)
lat0 = deg2rad(ll0(:,1));
lon0 = deg2rad(ll0(:,2));

lat1 = deg2rad(ll1(:,1));
lon1 = deg2rad(ll1(:,2));

hav = haversine(lat1-lat0) + cos(lat0).*cos(lat1).*haversine(lon1-lon0);

d = 2 .* geoC.WGS84_A .* asin(sqrt(hav));
end


