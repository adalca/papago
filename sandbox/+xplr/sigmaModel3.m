%% Description
% We compare the model3-sigma, especially comparing against model1-sigma. 
% This is all done on real data from ADNI or buckner.
%
% Initial conclusions are discussed in onenote on 9/14/2015.
% Updates have been made after exploring sigma reconstruction methods (see xplorSigmaReconstruction)
%
% TODO: sample from various sigmas.

%% setup
o3 = ones(1, 3);
patchSize = o3 * 9;
location = LOC_VENTRICLE_EDGE+10; 
locpad = o3 * 1;
minW = 0.000000001;

doload = false;
dosample = false;
model3sigmaregType = 1;

% load data
[patchescell, vols] = subspacetools.loadpatches('niistruct', 'data/bucknerNiis_5_5.mat', 'niis', {'iso', 'mask', 'ds'});
[isopatches, maskpatches, dspatches] = patchescell{:}; 
% bisopatches = isopatches .* maskpatches + 

%% loading data
% if doload
%     % extract patches
%     if ~exist('niis', 'var'); load('data/bucknerNiis_5_5.mat'); end
%     [maskpatches, maskpatchidx] = subspacetools.nii2patchcol(niis.mask, patchSize, location, locpad);
%     [isopatches, isopatchidx] = subspacetools.nii2patchcol(niis.iso, patchSize, location, locpad);
%     [dspatches, dspatchidx] = subspacetools.nii2patchcol(niis.ds, patchSize, location, locpad);
% 
%     blurvols = cellfunc(@(x) volblur(x.img, 2), niis.iso);
%     [bisopatches, bisopatchidx] = subspacetools.vols2patchcol(blurvols, patchSize, location, locpad);
% end

%% Isotropic analyses
X = isopatches;
W = max(maskpatches, minW);
xc = bsxfun(@minus, X, mean(X));

% standard covariance of all isotropic voxels
isostdcov = cov(xc);
[~, isostdcoverr] = cholcov(isostdcov);

% model 3 sigma with full weights.
[isom3sw1, isom3sw1err] = model3sigma(xc, W*0+1);

% force model 1 sigma
[isom1s, isom1sw1err] = model1sigma(xc, W);

% model 3 sigma if isotropic sigma
xc = bsxfun(@minus, X, wmean(X, W));
[isom3s, isom3serr, isom3wt, isom3r] = model3sigma(xc, W, model3sigmaregType);

% model 3 sigma with permutted W
rW = reshape(W(randperm(numel(W))), size(W));
xc = bsxfun(@minus, X, wmean(X, rW));
[isom3srW, isom3srWerr] = model3sigma(xc, rW);

isosigmas = {isostdcov, isom1s, isom3sw1, isom3s, isom3srW};
isotitles = {'isocov', 'iso model1', 'iso model3 W=1', 'iso model3 realW', 'iso model3 realW-randpermed'};

fact = prod(locpad*2+1)*5;
ww = min(isom3wt, fact) ./ fact;
isom3sr = ww .* isom3s + (1-ww) .* isom3r;
wtimg = ww ./ max(ww(:)) * prctile(isom3s(:), 99);
isosigmas2 = {isom3r, wtimg, isom3sr, abs(isom3s - isostdcov), abs(isom3sr - isostdcov)};
isotitles2 = {'iso model3 realW regularizer', 'iso model3 realW (renormed)wt', 'min(w,1) * m3rW + (1-min(w,1)) * m3rWreg', '', ''};

%% dspatches using blurred volumes (not real ds patches) and rand 0/1 weights
% set up 0/1 W
percentInformative = 0.2;
% minW = 0.01;
W = (unifrnd(0, 1, size(X)) < percentInformative)*(1-minW) + minW;
% W = unifrnd(0.0001, 1, size(X));
dspatches = W .* isopatches + (1 - W) .* bisopatches;
dsxc = bsxfun(@minus, dspatches, mean(dspatches));

% force model 1 sigma
[wds01m1s, wds01m1sw1err] = model1sigma(dsxc, W);

% model 3 sigma of fakew-ds data
[wds01m3s, wds01m1serr] = model3sigma(dsxc, W);

% covariance of fakew-ds
wds01stdcov = cov(dsxc);

wds01sigmas = {wds01m1s, wds01m3s, abs(wds01m3s-isostdcov), wds01stdcov, abs(wds01stdcov-isostdcov)};
wds01titles = {'ds W=B(01) model1', 'ds W=B(01) model3', '|model3-isocov|', 'ds W=B(01) cov', '|ds W=B(01) cov-isocov|'};
wds01titles = cellfunc(@(x) [x, sprintf(' (%%inf: %3.2f)', percentInformative*100)], wds01titles);

%% dspatches using blurred volumes (not real ds patches) and rand weights
% set up 0/1 W
% minW = 0.01;
W = unifrnd(minW, 1, size(X));
percentInformative = sum(W(:)) ./ numel(W);
dspatches = W .* isopatches + (1 - W) .* bisopatches;
dsxc = bsxfun(@minus, dspatches, mean(dspatches));

% force model 1 sigma
[wdsm1s, wdsm1sw1err] = model1sigma(dsxc, W);

% model 3 sigma of fakew-ds data
[wdsm3s, dsm1serr] = model3sigma(dsxc, W);

% covariance of fakew-ds
wdsstdcov = cov(dsxc);

wdssigmas = {wdsm1s, wdsm3s, abs(wdsm3s-isostdcov), wdsstdcov, abs(wdsstdcov-isostdcov)};
wdstitles = {'ds W=r(0..1) model1', 'ds W=r(0..1) model3', '|model3-isocov|', 'ds W=r(0..1) cov',  '|ds W=r(0..1) cov - isocov|'};
wdstitles = cellfunc(@(x) [x, sprintf(' (%%inf: %3.2f)', percentInformative*100)], wdstitles);

%% dspatches using blurred volumes (not real ds patches) and real weights
W = max(maskpatches, minW);
dspatches = W .* isopatches + (1 - W) .* bisopatches;
percentInformative = sum(W(:)) ./ numel(W);
dsxc = bsxfun(@minus, dspatches, mean(dspatches));

% model 3 sigma of ds data
[dsm3s, dsm3serr, dsm3w, dsm3r] = model3sigma(dsxc, W, model3sigmaregType);
ww = min(dsm3w, fact) ./ fact;
dsm3comb = ww .* dsm3s + (1-ww) .* dsm3r;

% covariance of ds
dsstdcov = cov(dsxc);
[~, dsstdcoverr] = cholcov(dsstdcov);

% combine m3 and dscov
ww = (W' * W);
% ww = ww ./ max(ww(:)); % or ./ prod(locpad) * 3 with min(..., 1);
ww = min(ww ./ (prod(locpad*2+1) * 3), 1);
ddf = dsm3s .* ww + dsstdcov .* (1 - ww);

blurdssigmas = {dsm3s, abs(dsm3s - isostdcov), dsstdcov, ddf, dsm3comb};
blurdstitles = {'ds m3', '|ds m3 - isocov|', 'ds cov', 'weighted combo of ds m3 + ds cov', 'combo of m3 + reg'};
blurdstitles = cellfunc(@(x) [x, sprintf(' (%%inf: %3.2f)', percentInformative*100)], blurdstitles);


%% analyze
sigmas = [isosigmas, isosigmas2, wds01sigmas, wdssigmas, blurdssigmas];
titles = [isotitles, isotitles2, wds01titles, wdstitles, blurdstitles];
sigmasc = cellfunc(@(x) x(:), sigmas);
mm = minmax(cat(1, sigmasc{:})');
mm = max(0, mm);
mm = min(mm, prctile(cat(1, sigmasc{:})', 99));

view2D(sigmas, 'caxis', mm, 'titles', titles, 'subgrid', [5, 5]);
colormap gray;

% view2D({isostdcov, dsm3s, dsstdcov, ddf}, 'caxis', mm);
% colormap gray;

%% sample from dsm3s vs dsstdcov
if dosample
    % this cell is a mess right now :)
    % minW = 0.1;
    W = max(maskpatches, minW);
    dspatches = W .* isopatches + (1 - W) .* bisopatches;
    percentInformative = sum(W(:)) ./ numel(W);
    dsxc = bsxfun(@minus, dspatches, mean(dspatches));

    % model 3 sigma of ds data
    % [dsm3s, dsm3serr] = model3sigma(dsxc, W);
    % z = 0; ws = linspace(minW, 1, 100);
    % while dsm3serr > 0
    %     z = z + 1;
    %     w = max(maskpatches, ws(z));
    %     [dsm3s, dsm3serr] = model3sigma(dsxc, w);
    %     dsm3s = dsm3s + eye(size(dsm3s, 1))*0.002;
    %     [~, dsm3serr] = cholcov(dsm3s);
    % end
    % if z > 0, warning('Had to increase minW to %3.2f for sigma to be PD', ws(z)); end
    % dsm3serr

    % try another fix
    [dsm3s, ~] = model3sigma(dsxc, W);
    dsm3s = nearestSPD(dsm3s);
    % [V, D] = eig(dsm3so);
    % [s, si] = sort(diag(D), 'descend');
    % s(s < 0) = min(s(s > 0));
    % dsm3so = sorteda * diag(s) * invsorteda;

    % cs = cumsum(s);
    % sorteda = V(:, si);
    % cs = cs ./ max(cs);
    % f = find(cs > 0.95, 1, 'first');
    % invsorteda = inv(sorteda);
    % oo = sorteda(:, 1:f) * diag(s(1:f)) * invsorteda(1:f, :);
    % q = mvnrnd(zeros(1, f)', oo, 10);


    % covariance of ds
    dsstdcov = cov(dsxc);
    [~, dsstdcoverr] = cholcov(dsstdcov);


    dsstdgmm = wgmm(dspatches, W, 1, wmean(dspatches, W), dsstdcov, 1, inv(dsstdcov));
    ds3gmm = wgmm(dspatches, W, 1, wmean(dspatches, W), dsm3s, 1, inv(dsm3s));

    N = 20;
    [dsstdX, dsstdidx] = dsstdgmm.sample(N);
    [ds3X, ds3idx] = ds3gmm.sample(N);

    subspacetools.comparePatchGMMData(dsstdX, dsstdidx(:), ds3X, ds3idx(:), patchSize, N);
    subplot(1,2,1); title('dsstdgmm sample, data'); subplot(1,2,2); title('ds3gmm sample data');
end
