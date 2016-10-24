%% a rough script to quilt volumes.
% need to do this significantly more elegantly :).

setup;
PATH = '/data/vision/polina/projects/patchSynthesis/data/stroke/subvols/wholevol/mar12_2016/volrecon_iterppca_K5_ppcaK30_Ds7Us5_may05_v2/*.mat';
% PATH = '/data/vision/polina/projects/patchSynthesis/data/stroke/subvols/wholevol/mar12_2016/volrecon_iterppca_K5_ppcaK1-3-15_Ds7Us5_may06/*.mat';
dspath = '/data/vision/polina/projects/stroke/work/patchSynthesis/data/stroke/subvols/wholevol/mar12_2016/stroke_wholevol_Ds7Us5Reg_volumes.mat';
maskpath = '/data/vision/polina/projects/stroke/work/patchSynthesis/data/stroke/subvols/wholevol/mar12_2016/stroke_wholevol_Ds7Us5RegMask_volumes.mat';

PATH = '/data/vision/polina/projects/patchSynthesis/data/ADNI_T1_baselines/subvols/wholevol/mar12_2016/volrecon_iterppca_K5_ppcaK1-5-16_ppcaIter25_may12/*.mat';
PATH = '/data/vision/polina/projects/patchSynthesis/data/ADNI_T1_baselines/subvols/wholevol/mar12_2016/volrecon_iterppca_K5_ppcaK5_ppcaiter40_may13/*.mat';
PATH = '/data/vision/polina/projects/patchSynthesis/data/ADNI_T1_baselines/subvols/wholevol/mar12_2016/volrecon_iterppca_K5_ppcaK5-30-35_ppcaiters100_may14/*.mat';
dspath = '/data/vision/polina/projects/stroke/work/patchSynthesis/data/ADNI_T1_baselines/subvols/wholevol/mar12_2016/ADNI_T1_baselines_wholevol_Ds5Us5Reg_volumes.mat';
maskpath = '/data/vision/polina/projects/stroke/work/patchSynthesis/data/ADNI_T1_baselines/subvols/wholevol/mar12_2016/ADNI_T1_baselines_wholevol_Ds5Us5RegMask_volumes.mat';



if ~exist('vid', 'var')
    vid = 3;
end

d = sys.fulldir(PATH);
fprintf('found %d subvolumes\n', numel(d));

tic
q = cell(1); 
for i = 1:numel(d), 
    q{i} = matfile(d(i).name); 
end
fprintf('done matfile-ing:%3.2f\n', toc);

tic
locs = cell(1); 
vols = cell(1); 
for i = 1:numel(d);
    locs{i} = q{i}.reconLoc + 2; 
    vols{i} = q{i}.reconVols(:,:,:,vid); 
    wts{i} = q{i}.reconWeight; 
end
fprintf('done gathering volumes:%3.2f\n', toc);

zv = ones(size(vols{1}));
[z, idx] = patchlib.vol2lib(zv, [9,9,9]);
sub = dimsplit(1, ind2subvec(size(zv), idx(:)))';
c = cellfunc(@(x) reshape(x, [9,9,9]), dimsplit(1, z)');
[~, wtfull] = patchlib.quiltIrregularPatches(sub, c);
wtfull = wtfull./max(wtfull(:));

wts = cell(1); 
for i = 1:numel(d), 
    wts{i} = wtfull;
end

% get volume
tic
vol = patchlib.quiltIrregularPatches(locs, vols, 'weightPatches', wts);
fprintf('done quilting:%3.2f\n', toc);

% compare
q = matfile(dspath);
dsvol = q.volumes(:,:,:,vid);


q = matfile(dspath);
dsvol = q.volumes(:,:,:,vid);
qm = matfile(maskpath);
maskvol = qm.volumes(:,:,:,vid);

sz = size(dsvol);
vol(sz(1), sz(2), sz(3)) = 0;
combovol = maskvol .* dsvol + (1-maskvol) .* vol;

view3Dopt(dsvol, maskvol, vol, combovol);
