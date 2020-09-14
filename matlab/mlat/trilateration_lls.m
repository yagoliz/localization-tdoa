function [position] = trilateration_lls(n,pos,dist,W)
    size_dist = size(dist);
    if size_dist(1) == 1
        dist = dist';
    end
    A = pos - repmat(mean(pos),n,1);
    b = 0.5*(sum(pos.^2,2) - sum(mean(pos.^2)) - dist.^2 + mean(dist.^2));
%     position = (A'*W*A)^-1*A'*W*b;
%     Use more stable function and recover the mean that was substracted before
    position = linsolve(A'*W*A, A'*W*b) + mean(pos)';
end

