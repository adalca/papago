% get entropy of ds subvolumes
enSubvol = dsSubvols*nan; 
for i = 1:size(enSubvol, 4)
    enSubvol(:,:,:,i) = entropyfilt(dsSubvols(:,:,:,i), getnhood(strel('sphere', 2))); 
end
croppedEnSubvols = cropVolume(enSubvol, [diffPad + 1, 1], [subvolSize - diffPad, nSubj]);

% get patches
enPatchCol = patchlib.vol2lib(croppedEnSubvols, [atlPatchSize, 1]);
enDs = enPatchCol(dataridx, :);

% adni
enFit = polyfit([3.5, 5.0], [0.5, 0.85], 1);
wtThrEnDs = within([0.01, 0.85], polyval(enFit, enDs));

% stroke
enFit = polyfit([5.0, 6.0], [0.5, 0.85], 1);
wtThrEnDs = within([0.01, 0.85], polyval(enFit, enDs));

w = wtPatchCol(dataridx, :);
wtEnDs = w > wtThrEnDs;
fprintf('mean wt: %3.2f\n', mean(wtEnDs(:)))

%% visualize
idx = 8;
prep = @(x) cellfunc(@(x) reshape(double(x), atlPatchSize), dimsplit(1, x));
view3Dopt(cellfunc(@(x) prep(x(idx, :)), {yds, Ws, w, wtEnDs}))

%% wgmm setup
dataTmplDsEn = struct('Y', yds, 'W', double(wtEnDs), 'K', K);
   
% example setup
name = 'LS_diag_en_low35_high50';
data.(name) = dataTmplDsEn;
params.(name) = [lsopts, 'replicates', 1, 'init', 'wgmm', 'initArgs', struct('wgmm', wgDsLsDiag)];


% example setup
name = 'LS_dsIdxModel3wfit_en_low35_high50';
data.(name) = dataTmplDsEn;
params.(name) = [lsopts, 'init', 'latentSubspace-model3', 'initArgs', struct('wgmm', wgDsLs, 'sigmaCorr', 'wfit', 'sigmaCorrThr', 10)];
