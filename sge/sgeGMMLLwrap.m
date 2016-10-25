
irange = round(linspace(1, 18000, 1000));
krange = [2, 3, 5, 10, 15, 25];
clear ll;

ll = nan(numel(irange), numel(krange));
for ii = randperm(numel(irange))
    i = irange(ii);
    subvolfile = sprintf('/data/vision/polina/scratch/adalca/patchSynthesis/data//adni/subvols_v2//subvol_%d.mat', i);
    
    for ki = 1:numel(krange)
        k = krange(ki);    
        gmmfile = sprintf('/data/vision/polina/scratch/adalca/patchSynthesis/data//adni/gmms//gmm_%d_model0_K%d.mat', i, k);
        llfile = sprintf('/data/vision/polina/scratch/adalca/patchSynthesis/data//adni/tmp//ll_%d_model0_K%d.mat', i, k);
        
%         try
            % ll(ii, ki) = sgeGMMLL(gmmfile, subvolfile, 'brainIso2Ds5Us5sizeReg', 'true', '[9,9,9]', '5000', llfile);
            q = sgeGMMreconhack(gmmfile, subvolfile, 'brainIso2Ds5Us5sizeReg', 'true', '[9,9,9]', '1000', llfile);
            errs(ii, ki) = mean(q);
%         catch err
%             err
%         end
    end
    ii
    nfos{ii} = load(subvolfile, 'nfo');
    save errs2 errs;
end


%% see image
m = nan([nfos{1}.nfo.params.atlVolSize, numel(krange)]);
for ii = 1:numel(irange);
    q = nfos{ii}.nfo.gridloc;
    q = num2cell(q);
    m(q{:}, :) = reshape(ll(ii, :), [1, 1, 1, numel(krange)]);
end
[~, maxm] = max(m, [], 4);
view3D(maxm)