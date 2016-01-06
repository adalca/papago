%% impute
% this is an early analysis on whether we can impute data from building a pca model at each location
% for each patchcol.
%
% 07.23.2015

%% setup
% md should be loaded from patchstartup.m
patchSize = ones(1, 3) * 7;
location = [56, 72, 92];
locpad = ones(1, 3) * 2;
locsearch = [7, 7, 1]; %ones(1, 3) * 3;

%% load all the volumes
vi = verboseIter(1:md.getNumSubjects);
while vi.hasNext()
    [s, i] = vi.next();
    % Extract iso volume 
    isoniis{i} = md.loadModality('rigidRegBrain', i);
    % Extract ds volume
    dsniis{i} = md.loadModality('brainDsUs5Reg', i);
    % get mask
    maskniis{i} = md.loadModality('brainDsUs5RegMask', i);
end
vi.close();

%% Trial 1: completely separate optimization.
prange = arrayfunc(@(x, d, p) (x - d: x + d)', location, locsearch, patchSize);
prange = ndgrid2cell(prange{:})';
for p = 1:prod(locsearch * 2 + 1)
    tic;
    
    % get the location, the range for patche starts, and the volume needed to get those patches
    plocation = cellfun(@(r) r(p), prange);
    locRange = arrayfunc(@(x, d) (x-d:x+d)', plocation, locpad);
    volRange = arrayfunc(@(x, d, p) (x - d: x + d + p - 1)', plocation, locpad, patchSize);

    % load data
    cpatches = cell(1, md.getNumSubjects);
    isopatches = cell(1, md.getNumSubjects);
    mpatches = cell(1, md.getNumSubjects);
    centerpatchidx = cell(1, md.getNumSubjects);
   
     vi = verboseIter(1:md.getNumSubjects);
    while vi.hasNext()
        [s, i] = vi.next();

        % Extract iso patch
        nii = isoniis{i};
        subvol = nii.img(volRange{:});
        isopatches{i} = patchlib.vol2lib(subvol, patchSize);

        % Extract patch
        nii = dsniis{i};
        subvol = nii.img(volRange{:});
        cpatches{i} = patchlib.vol2lib(subvol, patchSize);

        % get mask
        masknii = maskniis{i};
        subvol = masknii.img(volRange{:}) > 0.1;
        mpatches{i} = patchlib.vol2lib(subvol, patchSize);

        centerpatchidx{i} = false(size(mpatches{i}, 1), 1);
        centerpatchidx{i}((size(centerpatchidx{i}, 1)+1)/2) = true;
    end
    vi.close();
    
    
    patches = cat(1, cpatches{:});
    isopatches = cat(1, isopatches{:});
    maskpatches = cat(1, mpatches{:});
    meanpatch = mean(patches);

    centeridx = cat(1, centerpatchidx{:});

    % learn separate patches
    K = 7;
    nRep = 50;
    tmppatches = patches;
    tmppatches(~maskpatches) = nan;

    [recovpatches, meanX] = pcaimpute(tmppatches, K, nRep, isopatches);

    recoveredrightpatches{p} = recovpatches(centeridx, :);
    isorightpatches{p} = isopatches(centeridx, :);
    dsrightpatches{p} = patches(centeridx, :);
    maskrightpatches{p} = maskpatches(centeridx, :);
    toc;
end

%% quilt each (sub)volume
for v = 1:md.getNumSubjects
    
    vpatches = cellfunc(@(p) p(v, :), recoveredrightpatches);
    vpatches = cat(1, vpatches{:});
    
    rpatches = cellfunc(@(p) p(v, :), isorightpatches);
    rpatches = cat(1, rpatches{:});
    
    dpatches = cellfunc(@(p) p(v, :), dsrightpatches);
    dpatches = cat(1, dpatches{:});
    
    mpatches = cellfunc(@(p) p(v, :), maskrightpatches);
    mpatches = cat(1, mpatches{:});
    
    isovolumes{v} = patchlib.quilt(rpatches, locsearch * 2 + 1, patchSize);
    dsvolumes{v} = patchlib.quilt(dpatches, locsearch * 2 + 1, patchSize);
    recovolumes{v} = patchlib.quilt(vpatches, locsearch * 2 + 1, patchSize); 
    maskvolumes{v} = patchlib.quilt(mpatches, locsearch * 2 + 1, patchSize); 
end

v = 7;
view3Dopt(volumes{v}, maskvolumes{v}, dsvolumes{v}, recovolumes{v})