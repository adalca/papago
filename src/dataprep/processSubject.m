function processSubject(md, subjid, dsRate, intensityNorm, atlMods, regType, padAmount, preregmod)
% Like processmd, but for a single subject.
% TODO switch from md applyfun style functions to the specific function
%
% regType = 'BRAIN' or 'ATLAS'

    if ischar(md), load(md); end
    if ischar(subjid) && ~isempty(str2double(subjid)), subjid = str2double(subjid); end
    if ischar(intensityNorm), intensityNorm = str2double(intensityNorm); end
    if ischar(atlMods), load(atlMods); end
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

    %% Special Case: simulate the DsXUs2-iso to Iso-DsXUs2-Ds2Us2
    niiname = [tempname, '.nii.gz'];
    IsoDsUs2size = sprintf('Iso2Ds%dUs%dsize', dsRate, 2);
    downsampleNii(md.getModality(IsoDsUs2size, subjid), [1, 1, 2], niiname);
    IsoDsUs2sizeDs2Us2 = sprintf('Iso2Ds%dUs%dsize_Ds2Us2', dsRate, 2);
    IsoDsUs2sizeDs2Us2Mask = sprintf('Iso2Ds%dUs%dsize_Ds2Us2Mask', dsRate, 2);
    upsampleNii(niiname, md.getModality(IsoDsUs2sizeDs2Us2, subjid), ...
        md.getModality(IsoDsUs2sizeDs2Us2Mask, subjid), 'linear', 0, [1, 1, 2], true);
    delete(niiname);

    %% Perform registration via DsXUsX rigid registration
    % TODO: NN versions
    if ~exist('preregmod', 'var')
        atlfile = eval(sprintf('atlMods.%s_DS%d_US%d', upper(regType), dsRate, dsRate));
        preregmod = sprintf('Ds%dUs%dRegMat', dsRate, dsRate);
        DsUsReg = sprintf('Ds%dUs%dReg', dsRate, dsRate);
        DsUs = sprintf('Ds%dUs%d', dsRate, dsRate);
        md.register(DsUs, atlfile, 'rigid', 'monomodal', ...
            'saveModality', DsUsReg, 'savetformModality', preregmod, 'include', subjid);
    end
    
    usRatesSorted = sort(usRates, 'descend');
    for usRate = usRatesSorted
        % prepare atlas file (for applying warp)
        atlfile = eval(sprintf('atlMods.%s_DS%d_US%d', upper(regType), dsRate, usRate));

        % modality names for this dsRate and usRate
        IsoDsUssize = sprintf('Iso2Ds%dUs%dsize', dsRate, usRate);
        DsUs = sprintf('Ds%dUs%d', dsRate, usRate);
        DsUsMark = sprintf('Ds%dUs%dMask', dsRate, usRate);

        DsUsReg = sprintf('Ds%dUs%dReg', dsRate, usRate);
        DsUsRegMask = sprintf('Ds%dUs%dRegMask', dsRate, usRate);
        Iso2DsUssizeReg = sprintf('Iso2Ds%dUs%dsizeReg', dsRate, usRate);

        % apply "dsXusX" registration to modality and to mask
        md.register(DsUs, atlfile, 'rigid', 'monomodal', ...
            'saveModality', DsUsReg, 'loadtformModality', preregmod, 'include', subjid);
        md.register(DsUsMark, atlfile, 'rigid', 'monomodal', ...
            'saveModality', DsUsRegMask, 'loadtformModality', preregmod, 'include', subjid);
        md.register(IsoDsUssize, atlfile, 'rigid', 'monomodal', ...
            'saveModality', Iso2DsUssizeReg, 'loadtformModality', preregmod, 'include', subjid);
        
        % segmentations
        if doseg
            DsUsRegSeg = sprintf('Ds%dUs%dRegSeg', dsRate, dsRate); % use the original Ds5Us5!
            md.register(CroppedSeg, atlfile, 'rigid', 'multimodal', ...
                'saveModality', DsUsRegSeg, 'loadtformModality', preregmod, 'registeredVolumeInterp', 'nearest', 'include', subjid);        end
    end
    
    %% special case: warp of Iso-DsXUs2-Ds2Us2
    atlfile = eval(sprintf('atlMods.%s_DS%d_US%d', upper(regType), dsRate, 2));
    IsoDsUs2sizeDs2Us2 = sprintf('Iso2Ds%dUs%dsize_Ds2Us2', dsRate, 2);
    IsoDsUs2sizeDs2Us2Mask = sprintf('Iso2Ds%dUs%dsize_Ds2Us2Mask', dsRate, 2);
    IsoDsUs2sizeDs2Us2Reg = sprintf('Iso2Ds%dUs%dsize_Ds2Us2Reg', dsRate, 2);
    IsoDsUs2sizeDs2Us2MaskReg = sprintf('Iso2Ds%dUs%dsize_Ds2Us2MaskReg', dsRate, 2);
    md.register(IsoDsUs2sizeDs2Us2, atlfile, 'rigid', 'multimodal', ...
        'saveModality', IsoDsUs2sizeDs2Us2Reg, 'loadtformModality', preregmod, 'include', subjid);
    md.register(IsoDsUs2sizeDs2Us2Mask, atlfile, 'rigid', 'multimodal', ...
        'saveModality', IsoDsUs2sizeDs2Us2MaskReg, 'loadtformModality', preregmod, 'include', subjid);

    %% mdInterpmatWarp(md, dsmod, dsusmaskmod, rmod, atlvol, regmod)
%     % Transform medicalDataset sparse-slice volumes according to sparse-slice interpolant matrix, as
%     % opposed to warping dsXusX images
% 
%     if ismember('mdInterpmatWarp', steps); % unfinished.
% 
%         dsmod = 'Ds5Us5Mask';
%         dsusmaskmod = 'Ds5Us5InterpMat';
%         rmod = 'Ds5Us5Regwcor';
%         atlvol = BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_BRAIN_PROC_DS5_US5;
%         regmod = 'Ds5Us5InterpReg';
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
%         % q = md.loadVolume('Ds5Us5Reg', i);
%         % qn = q; q(isnan(warpedVol)) = nan;
%         % qiso = md.loadVolume('Iso2Ds5Us5sizeReg', i);
%         % qiso(isnan(warpedVol)) = nan;
%         % view3Dopt(warpedVol, q, qn, qiso);
%     end

    %% Perform registration via iso rigid registration
%     if ismember('isoreg', steps);
%         % register original without downsampling Note: in some sense, this is "true" rigid registration,
%         %   since we used the true data to perform the registration.
%         md.register('proc', atlMods.BUCKNER_ATLAS_BRAIN_PROC, 'rigid', 'multimodal', ...
%             'saveModality', 'rigidReg', 'savetformModality', 'rigidRegMat');
%         usRatesSorted = sort(usRates, 'descend');
%         for usRate = usRatesSorted
%             % prepare atlas file (for applying warp)
%             atlfile = eval(sprintf('atlMods.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', dsRate, usRate));
% 
%             regDsUs = sprintf('regDs%dUs%d', dsRate, usRate);
%             regDsUsMask = sprintf('regDs%dUs%dMask', dsRate, usRate);
%             RegMat = 'rigidRegMat';
% 
%             % apply "iso" (better?) registration to modality and to mask
%             md.register(DsUs, atlfile, 'rigid', ...
%                 'saveModality', regDsUs, 'loadtformModality', RegMat);
%             md.register(DsUsMark, atlfile, 'rigid', ...
%                 'saveModality', regDsUsMask, 'loadtformModality', RegMat);
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
%                 eval(sprintf('cropped%dnii = md.loadModality(''cropped%d'', i);', dsRate, dsRate));
%                 eval(sprintf('cropped%d = Cropped%dnii.img;', dsRate, dsRate));
%                 eval(sprintf('cropped%d_hdr = Cropped%dnii.hdr;', dsRate, dsRate));
%                 vars = [vars, sprintf('cropped%d', dsRate), sprintf('cropped%d_hdr', dsRate)];
% 
%                 %  isotropic to the size of dsXusY
%                 eval(sprintf('Iso2Ds%dUs%dsizenii = md.loadModality(''Iso2Ds%dUs%dsize'', i);', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('Iso2Ds%dUs%dsize = Iso2Ds%dUs%dsizenii.img;', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('Iso2Ds%dUs%dsize_hdr = Iso2Ds%dUs%dsizenii.hdr;', dsRate, usRate, dsRate, usRate));
%                 vars = [vars, sprintf('Iso2Ds%dUs%dsize', dsRate, usRate), sprintf('Iso2Ds%dUs%dsize_hdr', dsRate, usRate)];
% 
%                 eval(sprintf('bCropped%d_sigma2 = volblur(Cropped%d, 2);', dsRate, dsRate));
%                 vars = [vars, sprintf('bCropped%d_sigma2', dsRate)];
% 
%                 eval(sprintf('Ds%dUs%dMasknii = md.loadModality(''Ds%dUs%dMask'', i);', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('Ds%dUs%dMask = Ds%dUs%dMasknii.img;', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('Ds%dUs%dMask_hdr = Ds%dUs%dMasknii.hdr;', dsRate, usRate, dsRate, usRate));
%                 vars = [vars, sprintf('Ds%dUs%dMask', dsRate, usRate), sprintf('Ds%dUs%dMask_hdr', dsRate, usRate)]; 
% 
%                 % registered stuff
% %                 rigidRegnii = md.loadModality('rigidReg', i);
% %                 rigidReg = rigidRegnii.img;
% %                 rigidReg_hdr = rigidRegnii.hdr;
% %                 vars = [vars, 'rigidReg', 'rigidReg_hdr'];
% % 
% %                 brigidReg_sigma2 = volblur(rigidReg, 2);
% %                 vars = [vars, 'brigidReg_sigma2', 'rigidReg_hdr'];
% 
%                 eval(sprintf('Ds%dUs%dRegnii = md.loadModality(''Ds%dUs%dReg'', i);', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('Ds%dUs%dReg = Ds%dUs%dRegnii.img;', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('Ds%dUs%dReg_hdr = Ds%dUs%dRegnii.hdr;', dsRate, usRate, dsRate, usRate));
%                 vars = [vars, sprintf('Ds%dUs%dReg', dsRate, usRate), sprintf('Ds%dUs%dReg_hdr', dsRate, usRate)]; 
% 
%                 eval(sprintf('Ds%dUs%dRegMasknii = md.loadModality(''Ds%dUs%dRegMask'', i);', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('Ds%dUs%dRegMask = Ds%dUs%dRegMasknii.img;', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('Ds%dUs%dRegMask_hdr = Ds%dUs%dRegMasknii.hdr;', dsRate, usRate, dsRate, usRate));
%                 vars = [vars, sprintf('Ds%dUs%dRegMask', dsRate, usRate), sprintf('Ds%dUs%dRegMask_hdr', dsRate, usRate)];
% 
%                 eval(sprintf('Iso2Ds%dUs%dsizeRegnii = md.loadModality(''Iso2Ds%dUs%dsizeReg'', i);', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('Iso2Ds%dUs%dsizeReg = Iso2Ds%dUs%dsizeRegnii.img;', dsRate, usRate, dsRate, usRate));
%                 eval(sprintf('Iso2Ds%dUs%dsizeReg_hdr = Iso2Ds%dUs%dsizeRegnii.hdr;', dsRate, usRate, dsRate, usRate));
%                 vars = [vars, sprintf('Iso2Ds%dUs%dsizeReg', dsRate, usRate), sprintf('Iso2Ds%dUs%dsizeReg_hdr', dsRate, usRate)];
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

    