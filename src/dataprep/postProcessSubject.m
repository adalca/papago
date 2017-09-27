function postProcessSubject(md, subjid, dsRate, usRates)
% processing after ANTs Registrations are done.

    if ischar(md), load(md); end
    if ischar(subjid) && ~isempty(str2double(subjid)), subjid = str2double(subjid); end
	if ischar(dsRate), dsRate = str2double(dsRate); end
    if ischar(usRates), usRates = str2double(usRates); end
    
    % decide whether to do segmentations
    doseg = sys.isfile(md.getModality('seg', subjid));
    
    %% Perform registration via DsXUsX rigid registration
    antsfile = md.getModality(sprintf('Ds%dUs%dANTsAffine', dsRate, dsRate), subjid);
    preregmod = sprintf('Ds%dUs%dRegMat', dsRate, dsRate);
    preregfile = md.getModality(preregmod, subjid);
    tform = antsAffine2tformAffine3d(antsfile);
    tform = tform.invert;

    % unfortunately the antsAffine2tformAffine3d doesn't seem to
    % return a transform that exactly matches :(
    [optimizer, metric] = imregconfig('monomodal');
    subjnii = md.loadModality(sprintf('Ds%dUs%d', dsRate, dsRate), subjid);
    dstNii = md.loadModality(sprintf('Ds%dUs%dANTsReg', dsRate, dsRate), subjid);
    d = dstNii.hdr.dime.pixdim(2:4);
    zi = imwarp(subjnii.img, tform, 'OutputView', imref3d(size(dstNii.img), d(2), d(1), d(3)));
    bb = imregtform(zi, dstNii.img, 'similarity', optimizer, metric);
    tform.T = tform.T * bb.T;

    % save tform
    save(preregfile, 'tform');
    
    usRatesSorted = sort(usRates, 'descend');
    for usRate = usRatesSorted
        % prepare atlas file (for applying warp)
        antsNii = md.loadModality(sprintf('Ds%dUs%dANTsReg', dsRate, usRate), subjid);

        % modality names for this dsRate and usRate
        IsoDsUssize = sprintf('Iso2Ds%dUs%dsize', dsRate, usRate);
        DsUs = sprintf('Ds%dUs%d', dsRate, usRate);
        DsUsMark = sprintf('Ds%dUs%dMask', dsRate, usRate);

        DsUsReg = sprintf('Ds%dUs%dReg', dsRate, usRate);
        DsUsRegMask = sprintf('Ds%dUs%dRegMask', dsRate, usRate);
        Iso2DsUssizeReg = sprintf('Iso2Ds%dUs%dsizeReg', dsRate, usRate);

        % apply "dsXusX" registration to modality and to mask
        md.register(DsUs, antsNii, 'affine', 'monomodal', ...
            'saveModality', DsUsReg, 'loadtformModality', preregmod, 'include', subjid);
        md.register(DsUsMark, antsNii, 'affine', 'monomodal', ...
            'saveModality', DsUsRegMask, 'loadtformModality', preregmod, 'include', subjid);
        md.register(IsoDsUssize, antsNii, 'affine', 'monomodal', ...
            'saveModality', Iso2DsUssizeReg, 'loadtformModality', preregmod, 'include', subjid);
        
        % segmentations
        if doseg
            DsUsSeg = sprintf('Ds%dUs%dSeg', dsRate, usRate);
            DsUsRegSeg = sprintf('Ds%dUs%dRegSeg', dsRate, usRate); % use the original Ds5Us5!
            md.register(DsUsSeg, antsNii, 'affine', 'multimodal', ...
                'saveModality', DsUsRegSeg, 'loadtformModality', preregmod, 'registeredVolumeInterp', 'nearest', 'include', subjid);        
        end
    end
