function [F, grad] = optimminl1l2_ecef(X,doa,sensors,combinations,epsilon)
    global err_prev

    if nargin < 5
        epsilon = mean(X) * 1e-6;
    end

    d = ecef_distance(sensors, X);
    t = zeros(size(combinations,1),1);
    si = combinations(:,1);
    sj = combinations(:,2);
    
    t = d(si) - d(sj);
    
    err = doa - t;
    wei = 1./sqrt(1+(err).^2./2);
    F = sum(wei.*(err.^2));
    
    grad = zeros(1,3);
    if nargout > 1
%         grad(1) = -2 * sum(err .* ((-(sensors(si,1)-X(1))./d(si)) - (-(sensors(sj,1)-X(1))./d(sj))));
%         grad(2) = -2 * sum(err .* ((-(sensors(si,2)-X(2))./d(si)) - (-(sensors(sj,2)-X(2))./d(sj))));
%         grad(3) = -2 * sum(err .* ((-(sensors(si,3)-X(3))./d(si)) - (-(sensors(sj,3)-X(3))./d(sj))));
        grad(1) = (optimminl1l2_ecef(X+[epsilon,0,0],doa,sensors,combinations)-optimminl1l2_ecef(X-[epsilon,0,0],doa,sensors,combinations))/(2*epsilon);
        grad(2) = (optimminl1l2_ecef(X+[0,epsilon,0],doa,sensors,combinations)-optimminl1l2_ecef(X-[0,epsilon,0],doa,sensors,combinations))/(2*epsilon);
        grad(3) = (optimminl1l2_ecef(X+[0,0,epsilon],doa,sensors,combinations)-optimminl1l2_ecef(X-[0,0,epsilon],doa,sensors,combinations))/(2*epsilon);
        
        err_prev = err;
    end
end