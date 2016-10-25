ecmiter = 5;

threshold = 0.5;
nSampl = 5000;
clear sigmast meanst sigmas0 means0 sigmas means means_iso sigmas_iso
for k = 1:gmmK
    kidx = mi == k;
    
    % get isotropic data
    Xorig = bucknerIsoPatchCol(kidx, :); 
    X0orig = bucknerDsPatchCol(kidx, :); 
    W = bucknerDsMaskPatchCol(kidx, :);
    
    % sample the data.
    f = find(kidx);
    if numel(f) > nSampl
        r = randsample(numel(f), nSampl);
        X0orig = X0orig(r, :);
        W = W(r, :);
        Xorig = Xorig(r, :);
    end
    
    % patch-wise normalize data.
    X = bsxfun(@minus, Xorig, mean(Xorig, 2));
    X0 = bsxfun(@minus, X0orig, mean(X0orig, 2));
    
    % basic stats from iso and ds data
    sigmast(:,:,k) = cov(X); meanst(k, :) = mean(X);    
    sigmas0(:,:,k) = cov(X0); means0(k, :) = mean(X0);
    
    % set up downsampled data with nans in any region with W < thr
    X0(W < threshold) = nan;
    
%     for j = 1:100
        tic;
        X0 = X0orig;
        X0(W < threshold) = nan;
        [u,s,v] = svd(cov(X0orig)); 
        %Winit = Winit .* randn(size(Winit)) .* std(Winit(:)) * 2;
%         [m,s] = ecmninitx(X0, 'twostage');
%         [u,s,v] = svd(s); 
        ds = diag(s);
        explained = cumsum(ds)./sum(ds);
        f = 5;
        explained(f)
        vinit = (1 ./ 719) * sum(diag(s((f+1):729, (f+1):729)));
        Winit = u(:, 1:f) * (sqrt(s(1:f, 1:f)) - vinit * eye(f));
        opts = struct('TolFun', 0.001, 'TolX', 0.001, 'Display', 'off', 'MaxIter', 10);
        [COEFF, SCORE, LATENT, MU, V, S] = ppcax(X0, f, 'Options', opts, 'W0', Winit);
        fprintf('ppca cluster %d done in %5.3fs\n', k, toc);
%         nloglk_new(j) = ppcax_incomplete_nlogl(X0, MU, S.W, S.v);
%         sdd(j) = sum(abs(S.Recon(:) - Xorig(:)));
%         cla; plot(nloglk_new, sdd, '+'); drawnow; 
%         hold on; plot(nloglk_new(j-1), sdd(j-1), 'or');
%         hold on; plot(nloglk_new(j), sdd(j), 'ok');
%     end
   
    % cov(S.Recon)
    
    % W = u(:, 1:10) * sqrt((s(1:10, 1:10)) - (1 ./ 719) * sum(diag(s(11:729, 11:729))) *
    % eye(10)); % from tipping abd bishop eq (7) and (8)
    % c = W * W';
    
    srecon = bsxfun(@minus, S.Recon, mean(S.Recon, 2));
    meanppca(k,:) = mean(srecon);
    sigmappca(:,:,k) = cov(srecon);
    
%     imagesc([sigmas_iso(:,:,k), wgs{k}.sigma(:,:,1), sigmas_iso(:,:,k) - wgs{k}.sigma(:,:,1)])
%     drawnow;

    tic;
    X0 = bsxfun(@minus, X0orig, mean(X0orig, 2));
    X0(W < threshold) = nan;
    [meansppcaecm(k, :), sigmasppcaecm(:,:,k), Z] = ecmnmlex(srecon, 'twostage', ecmiter, 0.001, meanppca(k,:), sigmappca(:,:,k) + eye(size(X,2))*0.00001); 
    fprintf('cluster %d done in %5.3fs\n', k, toc);
    
    
    tic;
    X0 = bsxfun(@minus, X0orig, mean(X0orig, 2));
    c = cov(X0);
    [u,s,v] = svd(c);
%     vinit = (1 ./ 719) * sum(diag(s((f+1):729, (f+1):729)));
%     Winit = u(:, 1:f) * (sqrt(s(1:f, 1:f)) - vinit * eye(f));
    s(6:end, 6:end) = 0;
    c = u*s*v'+ eye(size(X0, 2)) * 0.00001;
    X0(W < threshold) = nan;
    m = nanmean(X0);
    %error('stop. this is bad m, c init');
%     [means(k, :), sigmas(:,:,k), Z] = ecmnmlex(X0, 'twostage', ecmiter, 0.001, m, c, 5); 
    %[means(k, :), sigmas(:,:,k), Z] = ecmnmlex(X0, 'twostage', ecmiter, 0.001);
    fprintf('cluster %d done in %5.3fs\n', k, toc);

end


wgtppca = wgmm(meanppca, sigmappca + repmat(eye(size(X,2))*0.00001, [1,1,gmmK]), hist(mi, 1:gmmK)./sum(mi));
[qselppca, reconLoc, cntvol] = papago.subvolRecon(wgtppca, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);


wgt = wgmm(means, sigmas + repmat(eye(size(X,2))*0.00001, [1,1,gmmK]), hist(mi, 1:gmmK)./sum(mi));
[qsel, reconLoc, cntvol] = papago.subvolRecon(wgt, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);


wgtppcaecm = wgmm(meansppcaecm, sigmasppcaecm + repmat(eye(size(X,2))*0.00001, [1,1,gmmK]), hist(mi, 1:gmmK)./sum(mi));
[qselppcaecm, reconLoc, cntvol] = papago.subvolRecon(wgtppcaecm, subvolLoc, subvolSize, atlPatchSize, crmethod, keepr, ...
    dsSubjInAtlNii.img, dsSubjInAtlMaskVol, dsSubjVol, dsSubjWeightVol, atlLoc2SubjSpace, extraReconArg);

