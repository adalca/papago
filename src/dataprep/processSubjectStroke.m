function processSubjectStroke(md, subjid, intensityNorm, atlMods, preregmod)
% Like processmd, but for a single subject.
% Assumes images have gone through brain extraction and wmmatch

    if ischar(md), load(md); end
    if ischar(subjid) && ~isempty(str2double(subjid)), subjid = str2double(subjid); end
    if ischar(intensityNorm), intensityNorm = str2double(intensityNorm); end
    if ischar(atlMods), load(atlMods); end
	
    % estimate dsrate.
    brainnii = md.loadModality('brain', subjid);
    dsRate = round(brainnii.hdr.dime.pixdim(4) ./ brainnii.hdr.dime.pixdim(2));
    assert(dsRate == (round(brainnii.hdr.dime.pixdim(4) ./ brainnii.hdr.dime.pixdim(3))));
    assert(isclean(dsRate));
    
    % us rates
    usRates = 1:dsRate;
    
    %% normalize intensity transform images to be between 0 to 1, and crop to a bounding box 
    md.normalize('brain', intensityNorm, 'procBrain', 'include', subjid);
        
    % get bounding box, crop procBrain and re-save to cropped brain.
    % [~, ~, bbrange, ~] = md.boundingBox('procBrain', 'procBrain', 'include', subjid);
    mdBoundingBoxRange(md, 'procBrain', 'procBrain', 'include', subjid);
    
    %% copy to downsample
    brainDs = sprintf('brainDs%d', dsRate);
    nii = md.loadModality('procBrain', subjid);
    md.saveModality(nii, brainDs, subjid);

    %% upsample
    
    % DsUs - upsample Ds volumes to usRate
    for usRate = usRates
        % modality names for this dsRate and usRate
        brainDsIso = sprintf('brainDs%dIso%d', dsRate, usRate); % this is isotropic in slices, but downsampled in z.
        brainDsUs = sprintf('brainDs%dUs%d', dsRate, usRate);
        brainDsUsMark = sprintf('brainDs%dUs%dMask', dsRate, usRate);
        brainDsUsNN = sprintf('brainDs%dUs%dNN', dsRate, usRate);
        brainDsUsNNMark = sprintf('brainDs%dUs%dNNMask', dsRate, usRate);

        % prepare functions to upsample
        dsfn = @(x, y) downsampleNii(x, [dsRate/usRate, dsRate/usRate, 1], y, false, 'nn');
        usfnlin = @(x, y, m) upsampleNii(x, y, m, 'linear', 0, [1, 1, usRate], true);
        usfnnn = @(x, y, m) upsampleNii(x, y, m, 'nearest', 0, [1, 1, usRate], true);

        % first, downsample to isotropic low-quality size
        md.applyfun(dsfn, {brainDs, brainDsIso}, 'include', subjid); % meant to be upsampled to an isotropic-resolution (but bad quality) after ds.
        md.applyfun(usfnlin, {brainDsIso, brainDsUs, brainDsUsMark}, 'include', subjid);
        md.applyfun(usfnnn, {brainDsIso, brainDsUsNN, brainDsUsNNMark}, 'include', subjid);
    end
    
    %% Perform affine registration via DsXUsX rigid registration
    % TODO: NN versions
    if ~exist('preregmod', 'var')
        atlfile = eval(sprintf('atlMods.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', dsRate, dsRate));
        preregmod = sprintf('brainDs%dUs%dRegMat', dsRate, dsRate);
        brainDsUsReg = sprintf('brainDs%dUs%dReg', dsRate, usRate);
        md.register(brainDsUs, atlfile, 'rigid', 'multimodal', ...
            'saveModality', brainDsUsReg, 'savetformModality', preregmod, 'include', subjid);
    end
    
    usRatesSorted = sort(usRates, 'descend');
    for usRate = usRatesSorted
        % prepare atlas file (for applying warp)
        atlfile = eval(sprintf('atlMods.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', dsRate, usRate));

        % modality names for this dsRate and usRate
        brainDsUs = sprintf('brainDs%dUs%d', dsRate, usRate);
        brainDsUsMark = sprintf('brainDs%dUs%dMask', dsRate, usRate);

        brainDsUsReg = sprintf('brainDs%dUs%dReg', dsRate, usRate);
        brainDsUsRegMask = sprintf('brainDs%dUs%dRegMask', dsRate, usRate);

        % apply "dsXusX" registration to modality and to mask
        md.register(brainDsUs, atlfile, 'rigid', 'multimodal', ...
            'saveModality', brainDsUsReg, 'loadtformModality', preregmod, 'include', subjid);
        md.register(brainDsUsMark, atlfile, 'rigid', 'multimodal', ...
            'saveModality', brainDsUsRegMask, 'loadtformModality', preregmod, 'include', subjid);
    end
    