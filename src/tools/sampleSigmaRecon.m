function [x, c, clocal, sigma6vec] = sampleSigmaRecon(N, varargin)
% SAMPLESIGMARECON sample sigma reconstruction model
%   xspl = sampleSigmaRecon(N, dataSigma, triplet) test sampling from p(x3|x1) by going through x2.
%   N is the number of samples to produce. dataSigma is the large covariance matrix of all the data
%   involved (not necessarily just x1, x2 and x3). triplet is a 1x3 vector giving the indices of
%   [x1, x2, x3] in the original data. That is, dataSigma(triplet(1), triplet(2)) is the covariance
%   of x1 with x2. The algorithm then samples p(x2|x1) followed by p(x3|x2). Output xspl is a three
%   columned matrix with the samples of x1, x2 and x3, respectively.
%
%   x = sampleSigmaRecon(N, dataSigma, pair) specify only the pair of indices for x1 and x3. We will
%   then sample through all possible x2s, N samples for each. Warning: we currently include x1 and
%   x3 in xk
%
%   x = sampleSigmaRecon(N, sigmavec) specify the "sigma vector" instead. the "sigma vector" is a
%   vector holding the relevant covariances for this sampling: sigma is a Mx6, where each row has:
%   [sigma_22, sigma_12, sigma_11, sigma_33, sigma_23, sigma_22] 
%
%   [x, covar, covarlocal, sigma] = sampleSigmaRecon(...) also outputs the covariance of the
%   samples (x' * x), local covariance (i.e. covaraince of the samples for each row of sigma above).
%   Also returns the "sigma vector" created (or passed in)
%
% Contact: {adalca,klbouman}@csail.mit.edu

    % input checking
    narginchk(2, 3);

    % prepare sigma
    sigmaCovarsCell = {[2, 2], [1, 2], [1, 1], [3, 3], [2, 3], [2, 2]};
    pts2sigma = @(ds, t) arrayfun(@(u) ds(t(sigmaCovarsCell{u}(1)), t(sigmaCovarsCell{u}(2))), 1:6);
    if nargin == 2
        sigma6vec = varargin{1};
    else
        dataSigma = varargin{1};
        
        if numel(varargin{2}) == 2 % pair given only
            pair = varargin{2};
            nFeats = size(dataSigma, 1);
            sigma6vec = zeros(nFeats, 6);
            for k = 1:nFeats
                sigma6vec(k, :) = pts2sigma(dataSigma, [pair(1), k, pair(2)]);
            end
        else
            assert(numel(varargin{2}) == 3, 'third argument must be pair or triplet of points');
            triplet = varargin{2};
            sigma6vec = pts2sigma(dataSigma, triplet);
        end
    end

    % sample
    xz = cell(1, size(sigma6vec, 1));
    clocal = cell(size(sigma6vec, 1), 1);
    
    for i = 1:size(sigma6vec, 1)

        x1 = normrnd(0, sqrt(sigma6vec(i, 3)), N, 1);
        x2 = normrnd(sigma6vec(i, 2)/sigma6vec(i, 3)*x1, sqrt(sigma6vec(i, 1)-sigma6vec(i, 2)^2/sigma6vec(i, 3)));
        x3 = normrnd(sigma6vec(i, 5)/sigma6vec(i, 6)*x2, sqrt(sigma6vec(i, 4)-sigma6vec(i, 5)^2/sigma6vec(i, 6)));
        
        xz{i} = [x1, x2, x3];
        if size(sigma6vec, 1) > 1
            clocal{i} = [x1, x2, x3]' * [x1, x2, x3] / size(x1, 1);
        end
    end
    x = cat(1, xz{:});

    % compute covariance
    c = x' * x / size(x, 1);
    c = real(c); 
    assert(isclean(c), 'c is not clean');
    