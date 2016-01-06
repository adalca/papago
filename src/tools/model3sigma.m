function [sigma, sigmacore, sigmarecon, sigmamerge] = model3sigma(X, W, reconmethod, combomethod, fact)
% a quick wrapper for sigma model3

    % prepare options
    methods = struct('core', 'model3', 'recon', reconmethod, 'merge', combomethod);
    opt = struct('mergeargs', fact, 'sigmareg', 0);

    % prepare inputs
    K = 1;
    mu = zeros(K, size(X, 2));
    gammank = ones(size(X, 1), K);

    % run full sigma. this includes spd correction
    [sigma, ~, sigmacore, sigmarecon, sigmamerge] = ...
        wgmm.sigmafull(mu, X, W, K, gammank, methods, opt);
end
