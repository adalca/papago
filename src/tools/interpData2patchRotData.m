function [interpDataPatches, availableData, availableDataIdx] = ...
    interpData2patchRotData(interpData, subvolSize, diffPad, atlPatchSize, subjids)
% prepare R matrices for each patch from the overall interpData for the subvolumes.

    % compute useful parameters
    nSubj = numel(subjids);
    nPatchesPerVol = prod(patchlib.gridsize(subvolSize, atlPatchSize));
    idxPatchCol = patchlib.grid(subvolSize, atlPatchSize); 
    volSize = subvolSize + diffPad * 2;

    %   transform linear index into the top-left subscript of each patch in the first volume
    %   in the subject space
    atlPatchSub = ind2subvec(subvolSize, idxPatchCol(:)); % 
    atlPatchSub = dimsplit(2, bsxfun(@plus, atlPatchSub, diffPad)); 
    
    % get linear indexes in atlas space for each patch
    atlPatchIdx = cell(1, nPatchesPerVol);
    for pi = 1:nPatchesPerVol
        a = cellfunc(@(c) c(pi, :), atlPatchSub);
        a = cellfunc(@(c,q) c:c+q-1, a, dimsplit(2, atlPatchSize));
        a = ndgrid2cell(a{:}); % ndgrid of this patch
        atlPatchIdx{pi} = sub2ind(volSize, a{:}); % ind of this patch
        atlPatchIdx{pi} = atlPatchIdx{pi}(:)';
    end
    atlPatchIdx = cat(1, atlPatchIdx{:});


    % prepare variables.
    %   Note: storing Rall and Gall as large nPatchesPerVol-by-nSubj cells slows down MATLAB
    %   in general *a lot*. I'm unclear why and haven't been able to get an answer from
    %   MatlabCentral/StackXchange. Everyone seems perplexed.
    %   https://www.mathworks.com/matlabcentral/answers/312344-small-matrix-multiplication-much-slower-in-r2016b-than-r2016a
    interpDataPatches.allR = [];
    interpDataPatches.lenR = zeros(nPatchesPerVol, nSubj);
    %interpDataPatches.allG = [];
    %interpDataPatches.lenG = zeros(nPatchesPerVol, nSubj);
    interpDataPatches.yorigAll = cell(nPatchesPerVol, nSubj);
    interpDataPatches.yrotmasksAll = cell(nPatchesPerVol, nSubj);
    interpDataPatches.ydsmasksAll = cell(nPatchesPerVol, nSubj);
    interpDataPatches.boundingBoxes = cell(nPatchesPerVol, nSubj);
    interpDataPatches.completed = false(nPatchesPerVol, nSubj);

     % go through available subvolumes
    vic = verboseIter(subjids, 2);
    while vic.hasNext()
        vi = vic.next(); tic

        % some subvolumes are 
        if isempty(interpData.atl2SubjInterpMat{vi})
            fprintf(2, 'Skipping %d\n', vi);
            continue;
        end

        % setup local storing for R
        localR = cell(nPatchesPerVol, 1);
        %localG = cell(nPatchesPerVol, 1);

        % extract the relevant data 
        volR = interpData.atl2SubjInterpMat{vi};
        assert(numel(find(volR(:))) > 0, 'Subvolume R is blank');
        %volG = interpData.subj2AtlInterpMat{vi};
        volyorig = interpData.subVols{vi};
        volymask = interpData.subVolMasks{vi};
        assert(all(volymask(:) == 1 | volymask(:) == 0));
        volymask = logical(volymask);

        % go through each patch of this data
        for pi = 1:nPatchesPerVol
            r = volR(:, atlPatchIdx(pi, :));
            %g = volG(atlPatchIdx(pi, :), :);

            % select which 'subject space' voxels to keep
            sel = sum(r, 2)' > 0; %used to be >0.5 % use to also have | sum(g, 1) > 0.5;
            r = r(sel, :);
            %g = g(:, sel);

            % make a yorig mask (TODO: might be better to save this in the interpDataFile?)
            %   via vol2subvolInterpmat
            z = false(size(volyorig));
            z(sel) = true;
            if sum(sel(:)) == 0
                continue;
            end

            [z, ~, ~, interpDataPatches.boundingBoxes{pi, vi}] = boundingBox(z);
            % normalize: rows have to add up to 1.
            localR{pi} = bsxfun(@rdivide, r, sum(r, 2));
            %localG{pi} = bsxfun(@rdivide, g, sum(g, 2));
            interpDataPatches.yorigAll{pi, vi} = volyorig(sel);
            interpDataPatches.ydsmasksAll{pi, vi} = volymask(sel);
            interpDataPatches.rWeightAll{pi, vi} = sum(r, 2)';
            interpDataPatches.ydsmasksAllFullVoxels{pi, vi} = full(volymask(sel) & sum(r, 2)' > 0.97);
            interpDataPatches.yrotmasksAll{pi, vi} = z;

            interpDataPatches.lenR(pi, vi) = size(localR{pi}, 1);
            %interpDataPatches.lenG(pi, vi) = size(localG{pi}, 2);

            % mark as completed volume
            interpDataPatches.completed(pi, vi) = true;
        end
        rvals = cat(1, localR{:});
        assert(isclean(find(rvals(:))), 'R is not clean');
        % assert(isclean(cat(2, localG{:})), 'G is not clean');
        warning('Skipping G Everythere in this file since not necessary and costly. FIXME');

        % put back into cell array
        interpDataPatches.allR = [interpDataPatches.allR; cat(1, localR{:})];
        % interpDataPatches.allG = [interpDataPatches.allG, cat(2, localG{:})];

        fprintf('Parsing R variables for subvol took %3.2f\n', toc);
    end

    % clean up interpData -- a partial reason for MATLAB slowdown

    % TODO: create an "available data" matrix of which patches are available, and what the
    % start/ends are for those matrices. It's very confusing to work otherwise.
    availableData = false(nPatchesPerVol, nSubj);
    availableData(interpDataPatches.completed) = true;
    availableDataIdx = find(availableData(:));
    % assert(all(interpDataPatches.lenR(availableDataIdx) > 0));
    %assert(all(interpDataPatches.lenG(availableDataIdx) > 0));
    