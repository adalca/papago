function processSubjectStroke(md, subjid, intensityNorm, atlMods, preregmod)
% Like processmd, but for a single subject.
% Assumes images have gone through  extraction and wmmatch

    if ischar(md), load(md); end
    if ischar(subjid) && ~isempty(str2double(subjid)), subjid = str2double(subjid); end
    if ischar(intensityNorm), intensityNorm = str2double(intensityNorm); end
    if ischar(atlMods), load(atlMods); end
	
    % estimate dsrate.
    nii = md.loadModality('orig', subjid);
    dsRate = round(nii.hdr.dime.pixdim(4) ./ nii.hdr.dime.pixdim(2));
    assert(dsRate == (round(nii.hdr.dime.pixdim(4) ./ nii.hdr.dime.pixdim(3))));
    dsRate = str2double(dsRate);
    
    % us rates
    usRates = 1:dsRate;
    
    %% normalize intensity transform images to be between 0 to 1, and crop to a bounding box 
    md.normalize('orig', intensityNorm, 'proc', 'include', subjid);
        
    % get bounding box, crop proc and re-save to cropped .
    % [~, ~, bbrange, ~] = md.boundingBox('proc', 'proc', 'include', subjid);
    mdBoundingBoxRange(md, 'proc', 'proc', 'include', subjid);
    
    %% copy to downsample
    Ds = sprintf('Ds%d', dsRate);
    nii = md.loadModality('proc', subjid);
    md.saveModality(nii, Ds, subjid);

    %% upsample
    
    % DsUs - upsample Ds volumes to usRate
    for usRate = usRates
        % modality names for this dsRate and usRate
        DsUs = sprintf('Ds%dUs%d', dsRate, usRate);
        DsUsMark = sprintf('Ds%dUs%dMask', dsRate, usRate);
        DsUsNN = sprintf('Ds%dUs%dNN', dsRate, usRate);
        DsUsNNMark = sprintf('Ds%dUs%dNNMask', dsRate, usRate);

        % prepare functions to upsample
        dsfn = @(x, y) downsampleNii(x, [dsRate/usRate, dsRate/usRate, 1], y, false, 'nn');
        usfnlin = @(x, y, m) upsampleNii(x, y, m, 'linear', 0, [1, 1, usRate], true);
        usfnnn = @(x, y, m) upsampleNii(x, y, m, 'nearest', 0, [1, 1, usRate], true);

        % first, downsample to isotropic low-quality size
        md.applyfun(dsfn, {Ds, DsIso}, 'include', subjid); % meant to be upsampled to an isotropic-resolution (but bad quality) after ds.
        md.applyfun(usfnlin, {DsIso, DsUs, DsUsMark}, 'include', subjid);
        md.applyfun(usfnnn, {DsIso, DsUsNN, DsUsNNMark}, 'include', subjid);
    end
    
    %% Perform affine registration via DsXUsX rigid registration
    % TODO: NN versions
    if ~exist('preregmod', 'var')
        atlfile = eval(sprintf('atlMods.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', dsRate, dsRate));
        preregmod = sprintf('Ds%dUs%dRegMat', dsRate, dsRate);
        DsUsReg = sprintf('Ds%dUs%dReg', dsRate, usRate);
        DsUs = sprintf('Ds%dUs%d', dsRate, dsRate);
        md.register(DsUs, atlfile, 'rigid', 'multimodal', ...
            'saveModality', DsUsReg, 'savetformModality', preregmod, 'include', subjid);
    end
    
    usRatesSorted = sort(usRates, 'descend');
    for usRate = usRatesSorted
        % prepare atlas file (for applying warp)
        atlfile = eval(sprintf('atlMods.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', dsRate, usRate));

        % modality names for this dsRate and usRate
        DsUs = sprintf('Ds%dUs%d', dsRate, usRate);
        DsUsMark = sprintf('Ds%dUs%dMask', dsRate, usRate);

        DsUsReg = sprintf('Ds%dUs%dReg', dsRate, usRate);
        DsUsRegMask = sprintf('Ds%dUs%dRegMask', dsRate, usRate);

        % apply "dsXusX" registration to modality and to mask
        md.register(DsUs, atlfile, 'rigid', 'multimodal', ...
            'saveModality', DsUsReg, 'loadtformModality', preregmod, 'include', subjid);
        md.register(DsUsMark, atlfile, 'rigid', 'multimodal', ...
            'saveModality', DsUsRegMask, 'loadtformModality', preregmod, 'include', subjid);
    end
    