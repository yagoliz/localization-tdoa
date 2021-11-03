classdef geoC
    properties(Constant)
        % Constants
        % WGS84 ellipsoid Earth parameters
        WGS84_A = 6378137.0;
        WGS84_F = 1.0/298.257223563;
        WGS84_B = geoC.WGS84_A * (1 - geoC.WGS84_F);
        WGS84_ECC_SQ = 1 - geoC.WGS84_B * geoC.WGS84_B / (geoC.WGS84_A * geoC.WGS84_A);
        WGS84_ECC = sqrt(geoC.WGS84_ECC_SQ);

        % Average radius for a spherical Earth
        SPHERICAL_R = 6371e3;

        % Some derived values
        wgs84_ep = sqrt((geoC.WGS84_A^2 - geoC.WGS84_B^2) / geoC.WGS84_B^2);
        wgs84_ep2_b = geoC.wgs84_ep^2 * geoC.WGS84_B;
        wgs84_e2_a = geoC.WGS84_ECC_SQ * geoC.WGS84_A;
    end
end