function preProcessSubject(md, subjid, dsRate, intensityNorm, padAmount)
% See processSubject. 
% this does not do registration

    if ischar(md), load(md); end
    if ischar(subjid) && ~isempty(str2double(subjid)), subjid = str2double(subjid); end
    if ischar(intensityNorm), intensityNorm = str2double(intensityNorm); end
	if ischar(dsRate), dsRate = str2double(dsRate); end
    if ischar(padAmount), padAmount = str2double(padAmount); end
    
    % us rates
    usRates = 1:dsRate;
    
    % decide whether to do segmentations
    doseg = sys.isfile(md.getModality('seg', subjid));
    
    %% normalize intensity transform images to be between 0 to 1, and crop to a bounding box 
    md.normalize('orig', intensityNorm, 'proc', 'include', subjid);
    padNii(md.getModality('proc', subjid), md.getModality('proc', subjid), [], padAmount);
        
    % get bounding box, crop proc and re-save to cropped .
    % [~, ~, bbrange, ~] = md.boundingBox('proc', 'proc', 'include', subjid);
    mdBoundingBoxRange(md, 'proc', 'proc', 'include', subjid);
    
    % do the same to segmentation
    if doseg
        [~,~,bbrange] = boundingBoxNii(md.getModality('orig', subjid));
        nii = md.loadModality('seg', subjid);
        vol = nii.img(bbrange{:});
        md.saveModality(makeNiiLike(vol, nii), 'procSeg', subjid);
    end    
    
    %% fix pixel dimensions since we assume [1, 1, 1]
    % make sure pix dimentions are close enough
    i = subjid;
    h = md.loadModality('proc', i);
    pixdim = h.hdr.dime.pixdim(2:4);
    rpixdim = round(pixdim);
    if ~all(pixdim == rpixdim)
        assert(all(isclose(rpixdim, pixdim, 0.001)));
        fprintf('%d diff: %e\n', i, mean(abs(pixdim - rpixdim)));
        h.hdr.dime.pixdim(2:4) = rpixdim;
        md.saveModality(h, 'proc', i);
    end
    
    %% downsample
    % Ds - downsample to dsRate
    Ds = sprintf('Ds%d', dsRate);
    md.applyfun(@(x, y) downsampleNii(x, [1, 1, dsRate], y), {'proc', Ds}, 'include', subjid);
    % DsSeg - downsample Segmentation to dsRate
    if doseg
        DsSeg = sprintf('Ds%dSeg', dsRate);
        md.applyfun(@(x, y) downsampleNii(x, [1, 1, dsRate], y), {'procSeg', DsSeg}, 'include', subjid);        
    end

    %% re-upsample
    % DsUs - upsample Ds volumes to usRate
    for usRate = usRates
        % modality names for this dsRate and usRate
        DsIso = sprintf('Ds%dIso%d', dsRate, usRate); % this is isotropic in slices, but downsampled in z.
        DsUs = sprintf('Ds%dUs%d', dsRate, usRate);
        DsUsMark = sprintf('Ds%dUs%dMask', dsRate, usRate);
        DsUsNN = sprintf('Ds%dUs%dNN', dsRate, usRate);

        % prepare functions to upsample and downsample
        dsfn = @(x, y) downsampleNii(x, [dsRate/usRate, dsRate/usRate, 1], y, false, 'nn');
        usfnlin = @(x, y, m) upsampleNii(x, y, m, 'linear', 0, [1, 1, usRate], true);
        usfnnn = @(x, y, m) upsampleNii(x, y, m, 'nearest', 0, [1, 1, usRate], true);

        % first, downsample to isotropic low-quality size
        md.applyfun(dsfn, {Ds, DsIso}, 'include', subjid); % meant to be upsampled to an isotropic-resolution (but bad quality) after ds.
        md.applyfun(usfnlin, {DsIso, DsUs, DsUsMark}, 'include', subjid);
        md.applyfun(usfnnn, {DsIso, DsUsNN, DsUsMark}, 'include', subjid);
    end
    
    % crop iso to match DsXUsX
    pni = @(x, y, m, v) padNii(x, y, m, size(nii2vol(v)) - size(nii2vol(x)), 'post');
    Cropped = sprintf('cropped%d', dsRate);
    CroppedMask = sprintf('cropped%dMask', dsRate);
    DsDs = sprintf('Ds%dUs%d', dsRate, dsRate);
    md.applyfun(pni, {'proc', Cropped, CroppedMask, DsDs}, 'include', subjid);
    % crop iso seg to match DsXUsX
    if doseg
        CroppedSeg = sprintf('cropped%dSeg', dsRate);
        CroppedMask = sprintf('cropped%dMask', dsRate);
        md.applyfun(pni, {'procSeg', CroppedSeg, CroppedMask, sprintf('Ds%dUs%d', dsRate, dsRate)}, 'include', subjid);
    end

    % resize iso to match DsXUsX size
    for usRate = usRates
        % downsample cropped 
        dsfn = @(x, y) downsampleNii(x, [dsRate/usRate, dsRate/usRate, dsRate/usRate], y, false, 'nn'); 
        IsoDsUssize = sprintf('Iso2Ds%dUs%dsize', dsRate, usRate);
        md.applyfun(dsfn, {Cropped, IsoDsUssize}, 'include', subjid); 
        
        if doseg
            DsUsSeg = sprintf('Ds%dUs%dSeg', dsRate, usRate);
            md.applyfun(dsfn, {'procSeg', DsUsSeg}, 'include', subjid);        
        end
    end

    % check equality if iso-in-plane from planes of iso-ds'ed-images
    i = subjid;
    dsSubjMod = sprintf('Ds%dUs%d', dsRate, dsRate);
    dsSubjModMaskMod = sprintf('Ds%dUs%dMask', dsRate, dsRate);
    isoSubjMod = sprintf('Iso2Ds%dUs%dsize', dsRate, dsRate);
    subjiso = md.loadVolume(isoSubjMod, i);
    subjds = md.loadVolume(dsSubjMod, i);
    subjdsmask = md.loadVolume(dsSubjModMaskMod, i);
    err = abs(subjds - subjiso);
    fprintf('mean error: %f\n', mean(err(subjdsmask(:) > 0)))
