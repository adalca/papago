function processSubject(md, subjid, dsRate, intensityNorm, atlMods, preregmod)
% Like processmd, but for a single subject.
% TODO switch from md applyfun style functions to the specific function

    if ischar(md), load(md); end
    if ischar(subjid) && ~isempty(str2double(subjid)), subjid = str2double(subjid); end
    if ischar(intensityNorm), intensityNorm = str2double(intensityNorm); end
    if ischar(atlMods), load(atlMods); end

    % us rates
    usRates = 1:dsRate;
    
    %% normalize intensity transform images to be between 0 to 1, and crop to a bounding box 
    md.normalize('brain', intensityNorm, 'procBrain', 'include', subjid);
        
    % get bounding box, crop procBrain and re-save to cropped brain.
    % [~, ~, bbrange, ~] = md.boundingBox('procBrain', 'procBrain', 'include', subjid);
    mdBoundingBoxRange(md, 'procBrain', 'procBrain', 'include', subjid);
    % do the same to segmentation
    [~,~,bbrange] = boundingBoxNii(md.getModality('brain', subjid));
    nii = md.loadModality('seg', subjid);
    vol = nii.img(bbrange{:});
    md.saveModality(makeNiiLike(vol, nii), 'procBrainSeg', subjid);
    
    
    %% fix pixel dimensions since we assume [1, 1, 1]
    % make sure pix dimentions are close enough
    i = subjid;
    h = md.loadModality('procBrain', i);
    pixdim = h.hdr.dime.pixdim(2:4);
    rpixdim = round(pixdim);
    if ~all(pixdim == rpixdim)
        assert(all(isclose(rpixdim, pixdim, 0.001)));
        fprintf('%d diff: %e\n', i, mean(abs(pixdim - rpixdim)));
        h.hdr.dime.pixdim(2:4) = rpixdim;
        md.saveModality(h, 'procBrain', i);
    end
    
    %% downsample
    % Ds - downsample to dsRate
    brainDs = sprintf('brainDs%d', dsRate);
    md.applyfun(@(x, y) downsampleNii(x, [1, 1, dsRate], y), {'procBrain', brainDs}, 'include', subjid);
    % DsSeg - downsample Segmentation to dsRate
    brainDsSeg = sprintf('brainDs%dSeg', dsRate);
    md.applyfun(@(x, y) downsampleNii(x, [1, 1, dsRate], y), {'procBrainSeg', brainDsSeg}, 'include', subjid);        

    %% re-upsample
    % DsUs - upsample Ds volumes to usRate
    for usRate = usRates
        % modality names for this dsRate and usRate
        brainDsIso = sprintf('brainDs%dIso%d', dsRate, usRate); % this is isotropic in slices, but downsampled in z.
        brainDsUs = sprintf('brainDs%dUs%d', dsRate, usRate);
        brainDsUsMark = sprintf('brainDs%dUs%dMask', dsRate, usRate);
        brainDsUsNN = sprintf('brainDs%dUs%dNN', dsRate, usRate);
        brainDsUsNNMark = sprintf('brainDs%dUs%dNNMask', dsRate, usRate);

        % prepare functions to upsample and downsample
        dsfn = @(x, y) downsampleNii(x, [dsRate/usRate, dsRate/usRate, 1], y, false, 'nn');
        usfnlin = @(x, y, m) upsampleNii(x, y, m, 'linear', 0, [1, 1, usRate], true);
        usfnnn = @(x, y, m) upsampleNii(x, y, m, 'nearest', 0, [1, 1, usRate], true);

        % first, downsample to isotropic low-quality size
        md.applyfun(dsfn, {brainDs, brainDsIso}, 'include', subjid); % meant to be upsampled to an isotropic-resolution (but bad quality) after ds.
        md.applyfun(usfnlin, {brainDsIso, brainDsUs, brainDsUsMark}, 'include', subjid);
        md.applyfun(usfnnn, {brainDsIso, brainDsUsNN, brainDsUsNNMark}, 'include', subjid);
    end
    
    % crop iso to match DsXUsX
    pni = @(x, y, m, v) padNii(x, y, m, size(nii2vol(v)) - size(nii2vol(x)), 'post');
    brainCropped = sprintf('brainCropped%d', dsRate);
    brainCroppedMask = sprintf('brainCropped%dMask', dsRate);
    brainDsDs = sprintf('brainDs%dUs%d', dsRate, dsRate);
    md.applyfun(pni, {'procBrain', brainCropped, brainCroppedMask, brainDsDs}, 'include', subjid);
    % crop iso seg to match DsXUsX
    brainCroppedSeg = sprintf('brainCropped%dSeg', dsRate);
    brainCroppedMask = sprintf('brainCropped%dMask', dsRate);
    md.applyfun(pni, {'procBrainSeg', brainCroppedSeg, brainCroppedMask, sprintf('brainDs%dUs%d', dsRate, dsRate)}, 'include', subjid);

    % resize iso to match DsXUsX size
    for usRate = usRates
        % downsample cropped brain
        dsfn = @(x, y) downsampleNii(x, [dsRate/usRate, dsRate/usRate, dsRate/usRate], y, false, 'nn'); 
        brainIsoDsUssize = sprintf('brainIso2Ds%dUs%dsize', dsRate, usRate);
        md.applyfun(dsfn, {brainCropped, brainIsoDsUssize}, 'include', subjid); 
    end

    % check equality if iso-in-plane from planes of iso-ds'ed-images
    i = subjid;
    dsSubjMod = sprintf('brainDs%dUs%d', dsRate, dsRate);
    dsSubjModMaskMod = sprintf('brainDs%dUs%dMask', dsRate, dsRate);
    isoSubjMod = sprintf('brainIso2Ds%dUs%dsize', dsRate, dsRate);
    subjiso = md.loadVolume(isoSubjMod, i);
    subjds = md.loadVolume(dsSubjMod, i);
    subjdsmask = md.loadVolume(dsSubjModMaskMod, i);
    err = abs(subjds - subjiso);
    fprintf('mean error: %f\n', mean(err(subjdsmask(:) > 0)))

    %% Special Case: simulate the DsXUs2-iso to Iso-DsXUs2-Ds2Us2
    niiname = [tempname, '.nii.gz'];
    brainIsoDsUs2size = sprintf('brainIso2Ds%dUs%dsize', dsRate, 2);
    downsampleNii(md.getModality(brainIsoDsUs2size, subjid), [1, 1, 2], niiname);
    brainIsoDsUs2sizeDs2Us2 = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2', dsRate, 2);
    brainIsoDsUs2sizeDs2Us2Mask = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2Mask', dsRate, 2);
    upsampleNii(niiname, md.getModality(brainIsoDsUs2sizeDs2Us2, subjid), ...
        md.getModality(brainIsoDsUs2sizeDs2Us2Mask, subjid), 'linear', 0, [1, 1, 2], true);
    delete(niiname);

    %% Perform registration via DsXUsX rigid registration
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
        brainIsoDsUssize = sprintf('brainIso2Ds%dUs%dsize', dsRate, usRate);
        brainDsUs = sprintf('brainDs%dUs%d', dsRate, usRate);
        brainDsUsMark = sprintf('brainDs%dUs%dMask', dsRate, usRate);

        brainDsUsReg = sprintf('brainDs%dUs%dReg', dsRate, usRate);
        brainDsUsRegMask = sprintf('brainDs%dUs%dRegMask', dsRate, usRate);
        brainIso2DsUssizeReg = sprintf('brainIso2Ds%dUs%dsizeReg', dsRate, usRate);

        % apply "dsXusX" registration to modality and to mask
        md.register(brainDsUs, atlfile, 'rigid', 'multimodal', ...
            'saveModality', brainDsUsReg, 'loadtformModality', preregmod, 'include', subjid);
        md.register(brainDsUsMark, atlfile, 'rigid', 'multimodal', ...
            'saveModality', brainDsUsRegMask, 'loadtformModality', preregmod, 'include', subjid);
        md.register(brainIsoDsUssize, atlfile, 'rigid', 'multimodal', ...
            'saveModality', brainIso2DsUssizeReg, 'loadtformModality', preregmod, 'include', subjid);
        
        % segmentations
        brainDsUsRegSeg = sprintf('brainIso2Ds%dUs%dsizeRegSeg', dsRate, dsRate); % use the original Ds5Us5!
        md.register(brainCroppedSeg, atlfile, 'rigid', 'multimodal', ...
            'saveModality', brainDsUsRegSeg, 'loadtformModality', preregmod, 'registeredVolumeInterp', 'nearest', 'include', subjid);
    end
    
    %% special case: warp of Iso-DsXUs2-Ds2Us2
    atlfile = eval(sprintf('atlMods.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', dsRate, 2));
    brainIsoDsUs2sizeDs2Us2 = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2', dsRate, 2);
    brainIsoDsUs2sizeDs2Us2Mask = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2Mask', dsRate, 2);
    brainIsoDsUs2sizeDs2Us2Reg = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2Reg', dsRate, 2);
    brainIsoDsUs2sizeDs2Us2MaskReg = sprintf('brainIso2Ds%dUs%dsize_Ds2Us2MaskReg', dsRate, 2);
    md.register(brainIsoDsUs2sizeDs2Us2, atlfile, 'rigid', 'multimodal', ...
        'saveModality', brainIsoDsUs2sizeDs2Us2Reg, 'loadtformModality', preregmod, 'include', subjid);
    md.register(brainIsoDsUs2sizeDs2Us2Mask, atlfile, 'rigid', 'multimodal', ...
        'saveModality', brainIsoDsUs2sizeDs2Us2MaskReg, 'loadtformModality', preregmod, 'include', subjid);

    %% mdInterpmatWarp(md, dsmod, dsusmaskmod, rmod, atlvol, regmod)
%     % Transform medicalDataset sparse-slice volumes according to sparse-slice interpolant matrix, as
%     % opposed to warping dsXusX images
% 
%     if ismember('mdInterpmatWarp', steps); % unfinished.
% 
%         dsmod = 'brainDs5Us5Mask';
%         dsusmaskmod = 'brainDs5Us5InterpMat';
%         rmod = 'brainDs5Us5Regwcor';
%         atlvol = BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_BRAIN_PROC_DS5_US5;
%         regmod = 'brainDs5Us5InterpReg';
% 
%         atlnii = loadNii(atlvol);
% 
%         i = subjid;
%         dsSubjNii = md.getModality(dsmod, i);
%         dsusSubjmasknii = md.getModality(dsusmaskmod, i);
%         interpSubjFile = md.getModality(rmod, i);
%         regOutFile = md.getModality(regmod, i);
% 
%         paffine.warpvol(dsSubjNii, dsusSubjmasknii, interpSubjFile, ...
%             atlnii, regOutFile);
% 
%         % visualize if return
%         % q = md.loadVolume('brainDs5Us5Reg', i);
%         % qn = q; q(isnan(warpedVol)) = nan;
%         % qiso = md.loadVolume('brainIso2Ds5Us5sizeReg', i);
%         % qiso(isnan(warpedVol)) = nan;
%         % view3Dopt(warpedVol, q, qn, qiso);
%     end

    %% Perform registration via iso rigid registration
%     if ismember('isoreg', steps);
%         % register original without downsampling Note: in some sense, this is "true" rigid registration,
%         %   since we used the true data to perform the registration.
%         md.register('procBrain', atlMods.BUCKNER_ATLAS_BRAIN_PROC, 'rigid', 'multimodal', ...
%             'saveModality', 'rigidRegBrain', 'savetformModality', 'rigidRegMatBrain');
%         usRatesSorted = sort(usRates, 'descend');
%         for usRate = usRatesSorted
%             % prepare atlas file (for applying warp)
%             atlfile = eval(sprintf('atlMods.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', dsRate, usRate));
% 
%             regBrainDsUs = sprintf('regBrainDs%dUs%d', dsRate, usRate);
%             regBrainDsUsMask = sprintf('regBrainDs%dUs%dMask', dsRate, usRate);
%             brainRegMat = 'rigidRegMatBrain';
% 
%             % apply "iso" (better?) registration to modality and to mask
%             md.register(brainDsUs, atlfile, 'rigid', ...
%                 'saveModality', regBrainDsUs, 'loadtformModality', brainRegMat);
%             md.register(brainDsUsMark, atlfile, 'rigid', ...
%                 'saveModality', regBrainDsUsMask, 'loadtformModality', brainRegMat);
%         end
%     end
    

    %% copy some specific niftis to matfile
%     % TODO - save all the modalities to the matfile? (simpler code) the problem is this can be huge.
%     if ismember('matfile', steps);
%         for usRate = usRates;
%             vi = verboseIter(1:md.getNumSubjects);
%             while vi.hasNext();
%                 i = vi.next();
%                 matfilename = md.getModality('matfile', i);
%                 vars = {};
% 
%                 % original data stuff.
%                 eval(sprintf('brainCropped%dnii = md.loadModality(''brainCropped%d'', i);', dsRate, dsRate));
%                 eval(sprintf('brainCropped%d = brainCropped%dnii.img;', dsRate, dsRate));
%                 eval(sprintf('brainCropped%d_hdr = brainCropped%dnii.hdr;', dsRate, dsRate));
%                 vars = [vars, sprintf('brainCropped%d', dsRate), sprintf('brainCropped%d_hdr', dsRate)];
% 
%                 % brain isotropic to the size of dsXusY
%                 eval(sprintf('brainIso2Ds%dUs%dsizenii = md.loadModality(''brainIso2Ds%dUs%dsize'', i);', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('brainIso2Ds%dUs%dsize = brainIso2Ds%dUs%dsizenii.img;', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('brainIso2Ds%dUs%dsize_hdr = brainIso2Ds%dUs%dsizenii.hdr;', dsRate, usRate, dsRate, usRate));
%                 vars = [vars, sprintf('brainIso2Ds%dUs%dsize', dsRate, usRate), sprintf('brainIso2Ds%dUs%dsize_hdr', dsRate, usRate)];
% 
%                 eval(sprintf('bbrainCropped%d_sigma2 = volblur(brainCropped%d, 2);', dsRate, dsRate));
%                 vars = [vars, sprintf('bbrainCropped%d_sigma2', dsRate)];
% 
%                 eval(sprintf('brainDs%dUs%dMasknii = md.loadModality(''brainDs%dUs%dMask'', i);', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('brainDs%dUs%dMask = brainDs%dUs%dMasknii.img;', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('brainDs%dUs%dMask_hdr = brainDs%dUs%dMasknii.hdr;', dsRate, usRate, dsRate, usRate));
%                 vars = [vars, sprintf('brainDs%dUs%dMask', dsRate, usRate), sprintf('brainDs%dUs%dMask_hdr', dsRate, usRate)]; 
% 
%                 % registered stuff
% %                 rigidRegBrainnii = md.loadModality('rigidRegBrain', i);
% %                 rigidRegBrain = rigidRegBrainnii.img;
% %                 rigidRegBrain_hdr = rigidRegBrainnii.hdr;
% %                 vars = [vars, 'rigidRegBrain', 'rigidRegBrain_hdr'];
% % 
% %                 brigidRegBrain_sigma2 = volblur(rigidRegBrain, 2);
% %                 vars = [vars, 'brigidRegBrain_sigma2', 'rigidRegBrain_hdr'];
% 
%                 eval(sprintf('brainDs%dUs%dRegnii = md.loadModality(''brainDs%dUs%dReg'', i);', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('brainDs%dUs%dReg = brainDs%dUs%dRegnii.img;', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('brainDs%dUs%dReg_hdr = brainDs%dUs%dRegnii.hdr;', dsRate, usRate, dsRate, usRate));
%                 vars = [vars, sprintf('brainDs%dUs%dReg', dsRate, usRate), sprintf('brainDs%dUs%dReg_hdr', dsRate, usRate)]; 
% 
%                 eval(sprintf('brainDs%dUs%dRegMasknii = md.loadModality(''brainDs%dUs%dRegMask'', i);', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('brainDs%dUs%dRegMask = brainDs%dUs%dRegMasknii.img;', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('brainDs%dUs%dRegMask_hdr = brainDs%dUs%dRegMasknii.hdr;', dsRate, usRate, dsRate, usRate));
%                 vars = [vars, sprintf('brainDs%dUs%dRegMask', dsRate, usRate), sprintf('brainDs%dUs%dRegMask_hdr', dsRate, usRate)];
% 
%                 eval(sprintf('brainIso2Ds%dUs%dsizeRegnii = md.loadModality(''brainIso2Ds%dUs%dsizeReg'', i);', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('brainIso2Ds%dUs%dsizeReg = brainIso2Ds%dUs%dsizeRegnii.img;', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('brainIso2Ds%dUs%dsizeReg_hdr = brainIso2Ds%dUs%dsizeRegnii.hdr;', dsRate, usRate, dsRate, usRate));
%                 vars = [vars, sprintf('brainIso2Ds%dUs%dsizeReg', dsRate, usRate), sprintf('brainIso2Ds%dUs%dsizeReg_hdr', dsRate, usRate)];
%                 
%                 if sys.isfile(matfilename)
%                     save(matfilename, vars{:}, '-append', '-v7.3');
%                 else
%                     save(matfilename, vars{:}, '-v7.3');
%                 end
%             end
%             vi.close();
%         end
%     end

    