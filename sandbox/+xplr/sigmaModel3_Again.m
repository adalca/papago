reconmethod = 'greedy1';
combomethod = 'wfull-mult';


%% Looking at sigmas.

[m3dataSigmaIso, m3isoreg, m3isoFull] = model3sigma(isopatchesCen, maskpatches, reconmethod, combomethod);
[m3dataSigmaBisow, m3bisowreg, m3bisowFull] = model3sigma(bisowpatchesCen, maskpatches, reconmethod, combomethod); 

% visualize
figuresc();
p95 = prctile([dataSigma(:); m3dataSigmaIso(:); bdataSigma(:); m3dataSigmaBisow(:)], 95);
subplot(231); imagesc(dataSigma, [0, p95]); colormap gray; title('covariance of iso data');
subplot(232); imagesc(bdataSigma, [0, p95]); colormap gray; title('covariance of blurry data');
subplot(233); imagesc(maskSigma); colormap gray; title('mask "covariance"');
subplot(234); imagesc(m3dataSigmaIso, [0, p95]); colormap gray; title('model-3 cov of iso data');
subplot(235); imagesc(m3dataSigmaBisow, [0, p95]); colormap gray; title('model-3 cov of ds (bisow) data');
subplot(236); imagesc(m3bisowreg, [0, p95]); colormap gray; title('model-3 reg for ds (bisow) data');


% visualize sigmas
figuresc();
p95 = prctile([dataSigma(:); m3dataSigmaIso(:); bdataSigma(:)], 95);
subplot(231); imagesc(dataSigma, [0, p95]); colormap gray; title('covariance of iso data');
subplot(232); imagesc(m3dataSigmaIso, [0, p95]); colormap gray; title('model-3 cov of iso data');
subplot(233); imagesc(bdataSigma, [0, p95]); colormap gray; title('covariance of blurry data');
subplot(234); imagesc(m3isoreg, [0, p95]); colormap gray; title('iso model-3 sigma reg');

fact = prod(locpad*2+1)*15;
ww = min(maskSigma, fact) ./ fact;
mult = median(m3dataSigmaIso(ww == 1) ./ m3isoreg(ww == 1));
isom3sr = ww .* m3dataSigmaIso + (1-ww) .* m3isoreg .* mult;

subplot(235); imagesc(ww); colormap gray; title('mask "covariance" rescaled');
subplot(236); imagesc(isom3sr, [0, p95]); colormap gray; title('iso model-3 sigma heuristic combo');


% visualize sigmas
figuresc();
p95 = prctile([dataSigma(:); m3dataSigmaBisow(:); bdataSigma(:)], 95);
subplot(231); imagesc(dataSigma, [0, p95]); colormap gray; title('covariance of iso data');
subplot(232); imagesc(m3dataSigmaBisow, [0, p95]); colormap gray; title('model-3 cov of ds bisow data');
subplot(233); imagesc(bdataSigma, [0, p95]); colormap gray; title('covariance of blurry data');
subplot(234); imagesc(m3bisowreg, [0, p95]); colormap gray; title('bisow model-3 sigma reg');
fact = prod(locpad*2+1)*15;
ww = min(maskSigma, fact) ./ fact;
mult = median(m3dataSigmaBisow(ww == 1) ./ m3bisowreg(ww == 1));
isom3sr = ww .* m3dataSigmaBisow + (1-ww) .* m3bisowreg .* mult;
subplot(235); imagesc(ww); colormap gray; title('mask "covariance" rescaled');
subplot(236); imagesc(isom3sr, [0, p95]); colormap gray; title('bisow model-3 sigma combo');


%% determine the fact manually
figuresc();
facts = linspace(eps, 10, 100);
qerr = zeros(1, numel(facts));
for f = 1:numel(facts) 
    fact = prod(locpad*2+1)*facts(f);
    ww = min(maskSigma, fact) ./ fact;
    mult = median(m3dataSigmaBisow(ww == 1) ./ m3bisowreg(ww == 1));
    isom3sr = ww .* m3dataSigmaBisow + (1-ww) .* m3bisowreg .* mult;
    qerr(f) = mse(isom3sr(:), dataSigma(:));
end
plot(facts, qerr, '.-');
