%% xplorEffectOfLinearInterpolationOnDSData
% Explore the effect of linear interpolation on the planes of data 
% that we trust. In this file, we use the Buckner Data

% volumes and parameters
subjNo = 4;
nThr = 100;
blurSigma = 1.5;

thresholds = linspace(0, 1, 20);
load('data/bucknerNiis_5_5');
m = niis.mask{subjNo}.img;
dsvolume = niis.ds{subjNo}.img;
isovolume = niis.iso{subjNo}.img;
bisovolume = (1-m) .* volblur(niis.iso{subjNo}.img, blurSigma) + m .* isovolume;

%% compute
% go over each threshold, adn compute the MSD of DS and ISO that are within those voxels whose
% weights are at least thr.
dsdst = zeros(1, nThr);
bisodst = zeros(1, nThr);
sm = zeros(1, nThr);
for i = 1:nThr 
    ll = thresholds(i); 

    dsdst(i) = msd(dsvolume(m(:) >= ll), isovolume(m(:) >= ll)) ./ mean(isovolume(m(:) >= ll));
    bisodst(i) = msd(bisovolume(m(:) >= ll), isovolume(m(:) >= ll)) ./ mean(isovolume(m(:) >= ll));
    sm(i) = sum(m(:) > ll);
end

%% visualize
figure(7); 

subplot(121);
plot(thresholds, dsdst, '.-'); hold on;
plot(thresholds, bisodst, '.-');
axislabels('threshold', 'normed MSD', 'MSD of voxels within weights >= thresh')
legend({'DS', sprintf('(1-w)*bluriso + w*iso mix with sigma=%f', blurSigma)});

subplot(122); 
plot(thresholds, sm);
axislabels('threshold', 'Count', 'Count of weight voxels >= thresh')
