%% Explore ADNI Data 
% Explore ADNI data in terms of aspects relevant to the subspace problem

%% visualize volumes and masks
usRate = dsRate;

% load volumes from adni md
dsregnii = md.loadModality(sprintf('regBrainDs%dUs%d', dsRate, usRate), reconSubj);
isoregnii = md.loadModality('rigidRegBrain', reconSubj);
maskregnii = md.loadModality(sprintf('regBrainDs%dUs%dMask', dsRate, usRate), reconSubj);

% visualize
view3Dopt(dsregnii.img, isoregnii.img, maskregnii.img, ...
    dsregnii.img .* maskregnii.img, isoregnii.img .* maskregnii.img);

% See the in-plane data a little closer
view3Dopt(dsregnii.img .* maskregnii.img, isoregnii.img .* maskregnii.img);
% Note: seem fairly similar so not bad

%% xplorEffectOfLinearInterpolationOnDSData
% do an analysis similar to xplorEffectOfLinearInterpolationOnDSData (or fix that one)

%% load patch data
location = [53, 53, 53];
volRange = arrayfunc(@(x, p) (x: x + p - 1)', location, patchSize);

[isopatches, layeridx, volidx] = colcmd('rigidRegBrain', location, locpad);
tic; isopatchesc = subspacetools.vols2patchcol(vols, patchSize, location, locpad); toc;
maskpatches = colcmd('regBrainDs5Us5Mask', location, locpad); 
dspatches = colcmd('regBrainDs5Us5', location, locpad);


dsregnii = md.loadModality(sprintf('regBrainDs%dUs%d', dsRate, usRate), reconSubj);
isoregnii = md.loadModality('rigidRegBrain', reconSubj);

%% visualize the resulting patches.
dsatlasPatch = reshape(dspatches(volidx == reconSubj & layeridx == 0, :), patchSize); 
dsatlasPatch_fromVol = dsregnii.img(volRange{:});
isoatlasPatch = reshape(isopatches(volidx == reconSubj & layeridx == 0, :), patchSize); 
isoatlasPatch_fromVol = isoregnii.img(volRange{:});
maskatlasPatch = reshape(maskpatches(volidx == reconSubj & layeridx == 0, :), patchSize);

view3Dopt(dsatlasPatch, isoatlasPatch, maskatlasPatch, dsatlasPatch_fromVol, isoatlasPatch_fromVol)

