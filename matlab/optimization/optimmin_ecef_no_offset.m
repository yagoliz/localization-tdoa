function [F, grad] = optimmin_ecef_no_offset(X,doa,sensors,combinations,epsilon,mu)

    if nargin < 5
        epsilon = mean(X(1:3)) * 1e-6;
        mu = mean(X(4:end))*1e-6;
    end

    N = size(sensors,1);

    d = ecef_distance(sensors, X(1:3));
    si = combinations(:,1);
    sj = combinations(:,2);
    offsets = X(4:end)';
    
    t = d(si) - d(sj) + (offsets(si) - offsets(sj));
    
    err = (doa - t).^2;
    F = sum(err);
    
    grad = zeros(1,3);
    if nargout > 1
        grad(1) = (optimmin_ecef(X+[epsilon,0,0,zeros(1,N)],doa,sensors,combinations)-optimmin_ecef(X-[epsilon,0,0,zeros(1,N)],doa,sensors,combinations))/(2*epsilon);
        grad(2) = (optimmin_ecef(X+[0,epsilon,0,zeros(1,N)],doa,sensors,combinations)-optimmin_ecef(X-[0,epsilon,0,zeros(1,N)],doa,sensors,combinations))/(2*epsilon);
        grad(3) = (optimmin_ecef(X+[0,0,epsilon,zeros(1,N)],doa,sensors,combinations)-optimmin_ecef(X-[0,0,epsilon,zeros(1,N)],doa,sensors,combinations))/(2*epsilon);
        for ii = 1:N
            ov = zeros(1,N);
            ov(ii) = mu;
            grad(3+ii) = (optimmin_ecef(X+[0,0,0,ov],doa,sensors,combinations)-optimmin_ecef(X-[0,0,0,ov],doa,sensors,combinations))/(2*mu);
        end
    end
end