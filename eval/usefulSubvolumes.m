function [selidx, grididx] = usefulSubvolumes(atlnii, atlsegnii, blurSigma, blurWindow, blurSigma2, blurWindow2, gridSpacing)
% decide on the useful subvolumes given an atlas and a gridSpacing.
%
% atlsegnii = loadNii([SYNTHESIS_DATA_PATH, 'buckner/atlases/wholevol/buckner61_seg_proc.nii.gz']);
% atlnii = loadNii([SYNTHESIS_DATA_PATH, 'buckner/atlases/wholevol/buckner61_proc.nii.gz']);
% usefulSubvolumes(atlnii, atlsegnii, 5, ones(1,3)*41, 5, ones(1,3)*17, ones(1,3)*6)

    patchSize = [1,1,1];
    patchOverlap = patchSize - gridSpacing;

    % get volumes
    seg = nii2vol(atlsegnii);
    atl = nii2vol(atlnii);

    % dilate brain
    dilatedSeg1 = volblur(double(seg>0), blurSigma, blurWindow) > 0;
    
    % take the AND of dilated brain with rough thr segmentation of whole brain
    newSeg = atl > 0.1 & dilatedSeg1;
    CC = bwconncomp(newSeg);
    [~, mi] = max(cellfun(@(x) numel(x), CC.PixelIdxList));
    newSeg = false(size(newSeg)); newSeg(CC.PixelIdxList{mi}) = true;
    
    % dilated new segmentation
    dilatedSeg = volblur(double(newSeg), blurSigma2, blurWindow2) > 0;
    
    % get all grid indexes
    grididx = patchlib.grid(size(seg), patchSize, patchOverlap);
    dots = zeros(size(seg));
    dots(grididx) = 1;
    
    % get grid indices within the segmentation
    selidx = find(dots(:) & dilatedSeg(:));
    nDotsInSeg = numel(selidx);
    
    % visualization 
%     view3Dopt(dilatedSeg1 + atl, double(atl > 0.1), dilatedSeg + atl + dots);
    view3Dopt(dilatedSeg + atl + dots);
    fprintf('Number of dots in seg: %d\n', nDotsInSeg);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% old code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    % % dilation instead of bluring. This is slow.
    % elementSize = ones(1, 3) * (serad+2-1);
    % n = size2ndgridvec(elementSize);
    % n = n - serad - 1;
    % elem = reshape(sqrt(sum(n.^2, 2)), elementSize);
    % se=strel(elem); %ones(50,50,50)); 
    % dilatedSeg = imdilate(seg>0,se); 
    % 
    % % visualize grid with arrows. Unforutnately, arrows need to be >=1 or so to see them, which gets
    % % confusing. The original idea was to have arrows of size 0, which woul hopfully show up as
    % % just circles.
    % sub = patchlib.grid(size(seg), patchSize, patchOverlap, 'sub');
    % 
    % f.u = sub{1}(:);
    % f.v = sub{2}(:);
    % f.w = sub{3}(:);
    % 
    % f.x = sub{1}(:)*1;
    % f.y = sub{2}(:)*1;
    % f.z = sub{3}(:)*1;
