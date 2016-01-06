function logp = logmvnpdf(x, mu, sigma, sigmainv)
% log of multivariate normal pdf as used in wgmm.
%   logp = logmvnpdf(x, mu, sigma) x is N-by-D, mu is K-by-D and sigma is D-by-D-by-K. returned logp
%   is N-by-K.
%
%   logp = logmvnpdf(x, mu, sigma, sigmainv) allows for the specification of the inverse of the
%   sigmas.
%
% part of wgmm project
%
% See Also: wgmm, logmvnpdf

    % parse inputs
    narginchk(3, 4);
    [N, D] = size(x);
    [K, mD] = size(mu);
    [sD1, sD2, sK] = size(sigma);
    assert(D == mD);
    assert(sD1 == sD2 && sD1 == D);
    assert(K == sK);
   
    % compute logmvnpdf one by one accross the K components
    logp = zeros(N, K);
    for k = 1:K
        if nargin == 3
            logp(:, k) = logmvnpdf(x, mu(k, :), sigma(:, :, k));
        else
            logp(:, k) = logmvnpdf(x, mu(k, :), sigma(:, :, k), sigmainv(:, :, k));
        end
    end
end
