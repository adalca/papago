%% (ssL) prepares the subvolume grid and the respective indices, given an atlas volume size.
% prepareSubvolumeIndices, uses usefulSubvolumes

outputPath = 'D:/Dropbox (MIT)/Research/patchSynthesis/data/buckner/subvols/wholevol/';
ds = 7;
us = 7;
gridSpacing = ones(1,3) * 10;

%% prep
gsText = sprintf('%dx%dx%d', gridSpacing);
gridSubFilename = [outputPath, sprintf('grididx2loc_ds%dus%d_gs%s.txt', ds, us, gsText)];
selSubFilename = [outputPath, sprintf('selidx2loc_ds%dus%d_gs%s.txt', ds, us, gsText)];

% get atlas files
atlsegnii = loadNii([SYNTHESIS_DATA_PATH, sprintf('buckner/atlases/wholevol/buckner61_seg_proc_ds%d_us%d.nii.gz', ds, us)]);
atlnii = loadNii([SYNTHESIS_DATA_PATH, sprintf('buckner/atlases/wholevol/buckner61_brain_proc_ds%d_us%d.nii.gz', ds, us)]);

%% compute locations

% parameters for "useful" subvolumes computation 
% (all subvolumes that touch the atlas brain a bit)
blurSigma = 5;
blurWindow1 = 41;
blurWindow2 = 21;
[selidx, grididx, brainmask] = usefulSubvolumes(atlnii, atlsegnii, ...
    blurSigma, ones(1, 3) * blurWindow1, blurSigma, ones(1, 3) * blurWindow2, gridSpacing);

% get grid locations
gridsub = ind2subvec(size(atlnii.img), grididx(:));

% get select grid locations
selsub = ind2subvec(size(atlnii.img), selidx(:));

%% save indices to files

% save to file
fid = fopen(gridSubFilename, 'w');
for i = 1:numel(grididx(:))
    fprintf(fid, '%d [%d,%d,%d]\n', grididx(i), gridsub(i, :));
end
fclose(fid);

% save to file
fid = fopen(selSubFilename, 'w');
for i = 1:numel(selidx(:))
    fprintf(fid, '%d [%d,%d,%d]\n', selidx(i), selsub(i, :));
end
fclose(fid);

%% closest to center
volLoc = [116, 104, 130]*us/ds;
% volLoc = [134, 197, 155]*us/ds;
% volLoc = [63, 42, 38];
nTop = 7^3;

% compute distances from volume location to each subvol
dists = sqrt(sum(bsxfun(@minus, volLoc, selsub).^2, 2));
[~, si] = sort(dists, 'ascend');

% split the top nTop and the rest into two files
[~, file, ~] = fileparts(selSubFilename);

fid = fopen(sprintf('%s/%s_top%d.txt', outputPath, file, nTop), 'w');
for i = 1:nTop
    fprintf(fid, '%d [%d,%d,%d]\n', selidx(si(i)), selsub(si(i), :));
end
fclose(fid);

fid = fopen(sprintf('%s/%s_rest%d.txt', outputPath, file, nTop), 'w');
for i = (nTop+1):numel(selidx(:))
    fprintf(fid, '%d [%d,%d,%d]\n', selidx(si(i)), selsub(si(i), :));
end
fclose(fid);

%% visualize with dots
brainidx = zeros(size(atlnii.img));
brainidx(selidx) = selidx;
braindots = double(brainidx > 0);

topidx = zeros(size(atlnii.img));
topidx(selidx(si(1:nTop))) = selidx(si(1:nTop));
topdots = double(topidx > 0);

view3Dopt(atlnii.img + topdots, topidx, atlnii.img + braindots, brainidx);
