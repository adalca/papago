% prepareSubvolumeIndices
PATH = 'D:/Dropbox (MIT)/Research/patchSynthesis/data/buckner/subvols/wholevol/';
gridSubFilename = [PATH, 'grididx2loc_ds7us5.txt'];
selSubFilename = [PATH, 'selidx2loc_ds7us5.txt'];
gridSpacing = ones(1,3)*7;
ds = 7;
us = 5;

% get
atlsegnii = loadNii([SYNTHESIS_DATA_PATH, 'buckner/atlases/wholevol/buckner61_seg_proc.nii.gz']);
atlnii = loadNii([SYNTHESIS_DATA_PATH, 'buckner/atlases/wholevol/buckner61_proc.nii.gz']);
atlsegnii = loadNii([SYNTHESIS_DATA_PATH, 'buckner/atlases/wholevol/buckner61_seg_proc_ds7_us5.nii.gz']);
atlnii = loadNii([SYNTHESIS_DATA_PATH, 'buckner/atlases/wholevol/buckner61_brain_proc_ds7_us5.nii.gz']);
[selidx, grididx] = usefulSubvolumes(atlnii, atlsegnii, 5, ones(1,3)*41, 5, ones(1,3)*17, gridSpacing);

%% save idx

% get grid locations
gridsub = ind2subvec(size(atlnii.img), grididx(:));
% save to file
fid = fopen(gridSubFilename, 'w');
for i = 1:numel(grididx(:))
    fprintf(fid, '%d [%d,%d,%d]\n', grididx(i), gridsub(i, :));
end
fclose(fid);

% get select grid locations
selsub = ind2subvec(size(atlnii.img), selidx(:));
% save to file
fid = fopen(selSubFilename, 'w');
for i = 1:numel(selidx(:))
    fprintf(fid, '%d [%d,%d,%d]\n', selidx(i), selsub(i, :));
end
fclose(fid);


%% closest to center
volLoc = [116, 104, 130]*us/ds;
nTop = 7^3;


df = sqrt(sum(bsxfun(@minus, volLoc, selsub).^2, 2));
[~, si] = sort(df, 'ascend');

% split the top nTop and the rest into two files
[~, file, ~] = fileparts(selSubFilename);

fid = fopen(sprintf('%s/%s_top%d.txt', PATH, file, nTop), 'w');
for i = 1:nTop
    fprintf(fid, '%d [%d,%d,%d]\n', selidx(si(i)), selsub(si(i), :));
end
fclose(fid);

fid = fopen(sprintf('%s/%s_rest%d.txt', PATH, file, nTop), 'w');
for i = (nTop+1):numel(selidx(:))
    fprintf(fid, '%d [%d,%d,%d]\n', selidx(si(i)), selsub(si(i), :));
end
fclose(fid);

%%
dots = zeros(size(atlnii.img));
dots(selidx(si(1:nTop))) = 1;
view3Dopt(atlnii.img + dots);

%%
dots = zeros(size(atlnii.img));
dots(selidx) = selidx;
view3Dopt(atlnii.img, dots);
