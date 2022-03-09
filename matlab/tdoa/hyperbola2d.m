function [ hyp ] = hyperbola2d(doa_meters, rx1, rx2, t)
% Function to generate the points of a hyperbola defined by the 2D position
% of 2 receivers and the tdoa value

    % Compute the hyperbola parameters
    c = norm(rx2-rx1)/2;
    
    % Check whether this is possible
    if abs(doa_meters)/2 > c
        disp(['<strong>TDOA delay (' num2str(doa_meters) ' meters) larger than RX distance (' num2str(c) ' meters) -> no solution possible </strong>']);
        doa_meters = sign(doa_meters) * 0.995 * c;
        disp(['<strong>ATTENTION: Correcting TODA delay to 0.995 * RX distance (maximum possible value) = ' num2str(0.995*c) '</strong>']);
    end    
    
    % Check again after correction
    if abs(doa_meters)/2 <= c
        % Compute midpoint and angle
        center = (rx2 + rx1)/2;
        theta = atan2((rx2(2)-rx1(2)),(rx2(1)-rx1(1)));

        R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

        % Computing the hyperbola parameters
        a = doa_meters/2;
        b = sqrt(c^2 - a^2);

        xpoints = a * cosh(t);
        ypoints = b * sinh(t);

        X_canonical = [fliplr(xpoints), xpoints; -fliplr(ypoints), ypoints];
        hyp = R * X_canonical + center';
    else
        disp('TDOA delay larger than RX distance -> no solution possible');
        hyp = [];
    end   
end