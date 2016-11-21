function [W, v, newSigma] = sigma2modelParams(sigma, dLow)
% extract subspace and residual variance from a full covariance

    % parse inputs
    narginchk(2, 2);
    dHigh = size(sigma, 1);
    assert(dHigh > dLow, 'subspace dimension is too large');

    % compute svd to get eigenspace    
    [u, s, ~] = svd(sigma);
    ds = diag(s);
    
    % lost variance
    v = 1 ./ (dHigh - dLow) * sum(ds(dLow+1:end));
    
    % new subspace
    W = u(:, 1:dLow) * sqrt(s(1:dLow, 1:dLow) - v * eye(dLow));
    
    % sigma based on new subspace
    newSigma = W * W' + v * eye(dHigh);
end
