function err = mccGMMreconhack(gmmfile, subvolfile, volnames, subtractpatchmean, patchSize, nSamples, randomseed, outfile)
% compute the reconstruction of a sampling of patches from isotropic data, destroyed at random, 
% in a subvolume under that subvolume's gmm
%
% MCC-friendly
%
%#ok<*ST2NM>
    
    % check inputs.
    narginchk(6, 8)
    [X, W] = subvol2XW(subvolfile, volnames, subtractpatchmean, patchSize, nSamples, randomseed);
    
    load(gmmfile, 'gmm');
    % get log-likehood of those patches
    % [~, ll] = wgmm.estep(X, W);
    err = zeros(1, nSamples);
    recond = cell(1, nSamples);
    for i = 1:nSamples
        x = X(i, W(i, :));
        mu = gmm.mu(:, W(i, :));
        sigma = gmm.sigma(W(i, :), W(i, :), :);
        pi = gmm.pi;
        gm = wgmm(mu, sigma, pi);
        gm.logpUpdateMethod = 'model0';
        p = gm.logpost(x, x*0+1);
        [~, mi] = max(p);
        xmu = X(i, :)' - gmm.mu(mi, :)';
        xmu(~W(i, :)) = nan; % set unknowns
        recond{i} = pcax.inpaint(xmu, gmm.sigma(:, :, mi)) + gmm.mu(mi, :)';
        err(i) = sqrt(sum((X(i, :) - recond{i}').^2));
    end
    caf;
    r = cat(2, recond{:})';
%     imagesc([r, X, abs(r-X)]); title(sprintf('mean(msd): %f mean(err): %f', mean(err), mean(abs(r(:)-X(:))))); colormap gray;
%     drawnow;
        

    % save/output
    %save(outfile, 'err');