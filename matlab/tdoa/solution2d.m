function [x, y] = solution2d(doa_array,sensors)
%SOLUTION2D Compute the TDOA algorithm and return the optimal position
%   Input:
%   - doa_array = array with the TDOAs between different pairs of sensors.
%   It is assumed that the first sensor is the reference and the rest are
%   TDOA calculations between the first and the corresponding sensor from 2
%   to M, being M the number of sensors
%   - sensors = matrix Mx2 with the positions of the sensors
%   
%   Output:
%   x, y = Estimation of the position of the unknown transmitter


% If we have 3 sensors, solution is exact (or there is none)
if length(sensors) == 3
    % This method is based on Fang's method:
    % https://ieeexplore.ieee.org/document/102710
    [x, y] = exactcomputation(doa_array, sensors);
   
    
% If we have 4 sensors, we can add a new variable and compute an exact
% solution
elseif length(sensors) == 4
    [A, b] = getMatrices(doa_array, sensors);
    result = A\b;
    x = result(2);
    y = result(3);
 

% Same as with 4 sensors, but we compute the least squares
elseif length(sensors) > 4
    [A, b] = getMatrices(doa_array, sensors);
    result = (A'*A)\A'*b;
    x = result(2);
    y = result(3);
    
else
    error('Not enough input sensors');
end

end

%% Get matrices
function [A, b] = getMatrices(doa_array, sensors)

% Let's construct the A and b matrices
A = zeros(length(doa_array),3);
b = zeros(length(doa_array),1);

for ii = 1:length(doa_array)
    % A matrix
    A(ii,1) = -doa_array(ii);
    A(ii,2) = sensors(1,1) - sensors(ii+1,1);
    A(ii,3) = sensors(1,2) - sensors(ii+1,2);
    
    % b array
    b(ii) = 0.5*(doa_array(ii)^2 + norm(sensors(1,:))^2 - norm(sensors(ii+1,:))^2);
    
end

end

%% Fang's method
function [x, y] = exactcomputation(doa_array, sensors)

% x = nan;
% y = nan;

% The core of this method is that reference is located in (0, 0)
% The second sensor is located on (b,0)
% For that we need to rotate the vectors
vectortoref = sensors(2:end,:) - sensors(1,:);

% With the basis line
basisline = vectortoref(1,:);
angle = atan2(basisline(2), basisline(1));
rotmat_0a = [cos(angle), -sin(angle); sin(angle), cos(angle)]; % Matrix to rotate from a to 0

% Now we put the sensors in the reference basis
sensorb = rotmat_0a' * basisline';
sensorc = rotmat_0a' * vectortoref(2,:)';

% We extract the coefficients
b = norm(sensorb);
c = norm(sensorc);
cx = sensorc(1);
cy = sensorc(2);

% We get the TDOA values
R12 = doa_array(1);
R13 = doa_array(2);

% Let's put some conditions to avoid errors
if abs(cy) > 1e-3
    if abs(R12) > 1e-3
        % We obtain the equation coefficients
        % y = g*x + h
        g = ((R13/R12)*b - cx)/cy;
        h = (c^2 - R13^2 + R12*R13*(1 - (b/R12)^2))/(2*cy);

        % 0 = d*x^2 + e*x + f
        d = -(1 + g^2 - (b/R12)^2);
        e = b*(1-(b/R12)^2) - 2*g*h;
        f = (R12^2/4)*(1-(b/R12)^2)^2 - h^2;

        % Now we obtain the positive term x and y
        xrotplus = (-e + sqrt(e^2 - 4*d*f))/(2*d);
        yrotplus = g*xrotplus + h;

        % Now we obtain the negative term x and y
        xrotminus = (-e - sqrt(e^2 - 4*d*f))/(2*d);
        yrotminus = g*xrotminus + h;
    
    else
        if abs(R13) > 1e-3
            x = b/2;
            alpha = (R13^2 - c^2 + 2*cx*x)/(2*R13);
            beta = cy/R13;
            
            d = beta^2 - 1;
            e = 2*alpha*beta;
            f = alpha^2 - x^2;
            
            % Now we obtain the positive term x and y
            xrotplus = x;
            yrotplus = (-e + sqrt(e^2 - 4*d*f))/(2*d);

            % Now we obtain the negative term x and y
            xrotminus = x;
            yrotminus = (-e - sqrt(e^2 - 4*d*f))/(2*d);

        else 
            % Add it to xrotminus (small trick to only output one value)
            xrotminus = b/2;
            yrotminus = (c^2 - cx*b)/(2*cy);
        end
    end

% end abs(cy) > 0
else
    if abs(cx*R12 - b*R13) < 1e-3
        warning('No solution possible');
        x = inf;
        y = inf;
        return
    else
        x = (R12^2*R13 - R12*R13^2+R12*cx^2 - R13*b^2) / ...
            (2*(cx*R12 - b*R13));
        
        xrotminus = x;
        xrotplus  = x;
        
        if abs(R12) > 1e-3
            yrotplus = sqrt((R12/2 - b^2/(2*R12) + (b*x)/(R12))^2 - x^2);
            yrotminus = -yrotplus;
            
        elseif abs(R12) < 1e-3 && abs(R13) > 1e-3
            yrotplus = sqrt((R13/2 - cx^2/(2*R13) + (cx*x)/(R13))^2 - x^2);
            yrotminus = -yrotplus;

        else
            warning('No solution possible');
            x = inf;
            y = inf;
            return
        end
        
    end
end

% Output the results in the proper basis
% Positive SQRT
if exist('xrotplus','var')
    resultplus = rotmat_0a * [xrotplus; yrotplus];

    % We add translate by the location of the sensor
    xplus = resultplus(1) + sensors(1,1);
    yplus = resultplus(2) + sensors(1,2);

    % We check the doa and see if the sign is equal to the doa_array
    plusvector = [real(xplus), real(yplus)];
    plusvectornorm = norm(plusvector - sensors(1,:)) - norm(plusvector - sensors(2,:));
    
    if abs(plusvectornorm) < 1e-3
        plusvectorsign = 0;
    else
        plusvectorsign = sign(plusvectornorm);
    end
    if plusvectorsign == sign(doa_array(1))
        x = xplus;
        y = yplus;
    end
end

% Negative SQRT
resultminus = rotmat_0a * [xrotminus; yrotminus];

% We add translate by the location of the sensor
xminus = resultminus(1) + sensors(1,1);
yminus = resultminus(2) + sensors(1,2);

% We check the doa and see if the sign is equal to the doa_array
minusvector = [real(xminus), real(yminus)];
minusvectornorm = norm(minusvector - sensors(1,:)) - norm(minusvector - sensors(2,:));

if abs(minusvectornorm) < 1e-3
    minusvectorsign = 0;
else
    minusvectorsign = sign(minusvectornorm);
end

if minusvectorsign == sign(doa_array(1))
    if exist('x','var')
        x = [x, xminus];
        y = [y, yminus];
    else
        x = xminus;
        y = yminus;
    end
end

if ~exist('x','var')
    x = 0;
    y = 0;
end

end
