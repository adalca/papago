function quiltedvolumes = volPatchImpute(volumes, validmasks, patchSize, inclpad, truevolumes)
% do function pcaimpute given a bunch of aligned (sub)volumes:
% build the library for each volume, and the sub-index for each location.
    nRep = 10;
    K = 7;

    
    
    [libs, ~, ~, gridsize] = patchlib.vol2lib(volumes, patchSize); % TODO: chancge truevolumes back to volumes
    [patchcol, mainidx] = libs2patchcells(libs, gridsize{1}, inclpad);
    
    [validlibs, ~, ~, gridsize] = patchlib.vol2lib(validmasks, patchSize);
    [validpatchcell, mainidx] = libs2patchcells(validlibs, gridsize{1}, inclpad);
    
    % S0. For each patch, get means, and fill in missing values (or maybe with linear interpolation instead?)
    % TODO: initlize with linear interpolation instead of means.
    for i = 1:numel(patchcol)
        missing = ~validpatchcell{i};
        z = patchcol{i};
        z(missing) = nan;
        m = repmat(nanmean(z), [size(z, 1), 1]);
        patchcol{i}(missing) = m(missing);
    end

    fprintf('%d, %f\n', ...
            0, mean(cellfun(@(v1, v2) msd(v1(:), v2(:)), volumes, truevolumes)));

    % Alternative optimization: Estimate PCA and Inpute
    for i = 1:nRep
        %   S1. For each patch, build PCA. Keep enough coefficients to cover 0.95% of variance?
        %       Note: in our tests, this might be bad, since you're probably just going to fill back in the same data.
        for c = 1:numel(patchcol)
            [coeff{c}, scores{c}, ~, ~, expl{c}, mu{c}] = pca(patchcol{c});
        end

%         % temporary: call pcaimpute.
%         [libs, ~, ~, gridsize] = patchlib.vol2lib(truevolumes, patchSize); % TODO: chancge truevolumes back to volumes
%         [realX, mainidx] = libs2patchcells(libs, gridsize{1}, inclpad);
%         q = patchcell{1}; q(~validpatchcell{1}) = nan;
%         pcaimpute(q, K, nRep, realX{1});
        
        % S2. For each patch, do either M1 or M2, then re-build libraries for S1.
        % M1: fill in data via PCA reconstruction, then quilt volumes from patches with filled in data
        for c = 1:numel(patchcol)
            missing = ~validpatchcell{c};
            Xi2 = pcax.recon(coeff{c}(:, 1:K), scores{c}(:, 1:K)', mu{c})';
            patchcol{c}(missing) = Xi2(missing);
        end
        % M2: Daniel Zoran's closed form solution?
        % TODO
        
        % go back from patch cells 2 libs, and index into libs
        selpatches = cellfunc(@(x, i) x(i, :), patchcol, mainidx);
        filledlibs = permutecellinception(selpatches);
        
        % quilt each volume
        quiltedvolumes = cellfunc(@(x) patchlib.quilt(x, gridsize{1}, patchSize), filledlibs);
        
        % build libraries. 
        [libs, ~, ~, gridsize] = patchlib.vol2lib(quiltedvolumes, patchSize);
        [patchcol, mainidx] = libs2patchcells(libs, gridsize{1}, inclpad);
        
        % TODO diplay some info?
        fprintf('%d, %f %f\n', ...
            i, mean(cellfun(@(v1, v2) msd(v1(:), v2(:)), quiltedvolumes, truevolumes')), ...
            mean(cellfun(@(x) sum(x(1:K)), expl)));
        disp('a');
    end
end


function [patchcell, mainidx] = libs2patchcells(libs, gridsize, inclpad)
% see subspace.nii2patchcol
    error('see subspace.nii2patchcol');
    
    patchcell = cell(prod(gridsize), 1);
    mainidx = cell(prod(gridsize), 1);

    % get "patch" library
    for i = 1:prod(gridsize)
        sub = ind2subvec(gridsize, i);
        range = arrayfunc(@(x, p, s) max(x - p, 1):min(x + p, s), sub, inclpad, gridsize);
        rangex = ndgrid2cell(range{:});
        idx = sub2indfast(gridsize, rangex{:});
        
        
        
        patches = cell(numel(libs), 1);
        indexes = cell(numel(libs), 1);
        for j = 1:numel(libs)
            patches{j} = libs{j}(idx(:), :);
            indexes{j} = idx(:) == i;
        end
        patchcell{i} = cat(1, patches{:});
        mainidx{i} = cat(1, indexes{:});
    end

end


function finalv = permutecellinception(c)
% go from a cell of prod(gridSize) elements, each with a [nVols x prod(patchSize)]  to a 
% cell of nVols elements, each with [prod(gridSize) x prod(patchSize)] elements
    
    v = cell(size(c{1}, 1), numel(c));
    for i = 1:numel(c)
        v(:, i) = mat2cell(c{i}, ones(size(c{i}, 1), 1), size(c{i}, 2));
    end
    
    finalv = cell(size(v, 1), 1);
    for i = 1:size(v, 1)
        finalv{i} =  cat(1, v{i, :});
    end
        
end
