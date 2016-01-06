function ll = mccGMMLL(gmmfile, subvolfile, volnames, subtractpatchmean, patchSize, nSamples, randomseed, outfile)
% compute the log likelihood of a sampling of patches in a subvolume under that subvolume's gmm
%
% MCC-friendly
%
%#ok<*ST2NM>
    
    % check inputs.
    narginchk(6, 8)
    [X, W] = subvol2XW(subvolfile, volnames, subtractpatchmean, patchSize, nSamples, randomseed);
    
    % get log-likehood of those patches via learned gmm
    load(gmmfile, 'gmm');
    ll = wgmm.logp(X, W);

    % save/output
    save(outfile, 'll');
    