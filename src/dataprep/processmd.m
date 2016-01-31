function processmd(md, dsRate, intensityNorm, atlmods, steps)
% PROCESSMD Process md dataset for a given downsampling rate for the restoration project
%
% the modalities should match the nadming from resotrationmd.
%
% dsRate = 5;
% intensityNorm = 255;
%
% TODO: have a parameter for forcing overwrite or skipping wherever
% possible/necessary, etc

    narginchk(4, 5);
    if nargin == 4
        steps = {'normalize', 'DsUs', 'isoreg', 'dsusreg', 'matfile', 'visualize'};
    end
    
    usRates = 1:dsRate;

    %% normalize intensity and size
    if ismember('normalize', steps);
        % transform images to be between 0 to 1, and crop to a bounding box 
        md.normalize('brain', intensityNorm, 'procBrain');
        
        [~, ~, bbrange, ~] = md.boundingBox('procBrain', 'procBrain');
    end

    %% downsample and upsample 
    if ismember('DsUs', steps);
        brainDs = sprintf('brainDs%d', dsRate);
        md.applyfun(@(x, y) downsampleNii(x, [1, 1, dsRate], y), {'procBrain', brainDs});

        for usRate = usRates
            % modality names for this dsRate and usRate
            brainDsIso = sprintf('brainDs%dIso%d', dsRate, usRate);
            brainDsUs = sprintf('brainDs%dUs%d', dsRate, usRate);
            brainDsUsMark = sprintf('brainDs%dUs%dMask', dsRate, usRate);
            brainDsUsNN = sprintf('brainDs%dUs%dNN', dsRate, usRate);
            brainDsUsNNMark = sprintf('brainDs%dUs%dNNMask', dsRate, usRate);

            % prepare functions to upsample and downsample
            dsfn = @(x, y) downsampleNii(x, [dsRate/usRate, dsRate/usRate, 1], y, false, 'nn');
            usfnlin = @(x, y, m) upsampleNii(x, y, m, 'linear', 0, [1, 1, usRate], true);
            usfnnn = @(x, y, m) upsampleNii(x, y, m, 'nearest', 0, [1, 1, usRate], true);

            % first, downsample to isotropic low-quality size
            md.applyfun(dsfn, {brainDs, brainDsIso}); % meant to be upsampled to an isotropic-resolution (but bad quality) after ds.
            md.applyfun(usfnlin, {brainDsIso, brainDsUs, brainDsUsMark});
            md.applyfun(usfnnn, {brainDsIso, brainDsUsNN, brainDsUsNNMark});
        end

        % crop iso to match DsXUsX
        pni = @(x, y, m, v) padNii(x, y, m, size(nii2vol(v)) - size(nii2vol(x)), 'post');
        brainCropped = sprintf('brainCropped%d', dsRate);
        brainCroppedMask = sprintf('brainCropped%dMask', dsRate);
        md.applyfun(pni, {'procBrain', brainCropped, brainCroppedMask, sprintf('brainDs%dUs%d', dsRate, dsRate)});

        % resize iso to match DsXUsX size
        for usRate = usRates
            % downsample cropped brain
            dsfn = @(x, y) downsampleNii(x, [dsRate/usRate, dsRate/usRate, dsRate/usRate], y, false, 'nn'); 
            brainIsoDsUssize = sprintf('brainIso2Ds%dUs%dsize', dsRate, usRate);
            md.applyfun(dsfn, {brainCropped, brainIsoDsUssize}); 
        end
    end

    %% Perform registration via iso rigid registration
    if ismember('isoreg', steps);
        % register original without downsampling Note: in some sense, this is "true" rigid registration,
        %   since we used the true data to perform the registration.
        md.register('procBrain', atlmods.BUCKNER_ATLAS_BRAIN_PROC, 'rigid', 'multimodal', ...
            'saveModality', 'rigidRegBrain', 'savetformModality', 'rigidRegMatBrain');
        usRatesSorted = sort(usRates, 'descend');
        for usRate = usRatesSorted
            % prepare atlas file (for applying warp)
            atlfile = eval(sprintf('atlmods.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', dsRate, usRate));

            regBrainDsUs = sprintf('regBrainDs%dUs%d', dsRate, usRate);
            regBrainDsUsMask = sprintf('regBrainDs%dUs%dMask', dsRate, usRate);
            brainRegMat = 'rigidRegMatBrain';

            % apply "iso" (better?) registration to modality and to mask
            md.register(brainDsUs, atlfile, 'rigid', ...
                'saveModality', regBrainDsUs, 'loadtformModality', brainRegMat);
            md.register(brainDsUsMark, atlfile, 'rigid', ...
                'saveModality', regBrainDsUsMask, 'loadtformModality', brainRegMat);
        end
    end
    
    %% Perform registration via DsXUsX rigid registration
    % TODO: NN versions
    if ismember('dsusreg', steps);
        usRatesSorted = sort(usRates, 'descend');
        for usRate = usRatesSorted
            % prepare atlas file (for applying warp)
            atlfile = eval(sprintf('atlmods.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', dsRate, usRate));

            % modality names for this dsRate and usRate
            brainIsoDsUssize = sprintf('brainIso2Ds%dUs%dsize', dsRate, usRate);
            brainDsUs = sprintf('brainDs%dUs%d', dsRate, usRate);
            brainDsUsMark = sprintf('brainDs%dUs%dMask', dsRate, usRate);
            %brainDsUsNN = sprintf('brainDs%dUs%dNN', dsRate, usRate);
            %brainDsUsNNMark = sprintf('brainDs%dUs%dNNMask', dsRate, usRate);

            brainDsUsReg = sprintf('brainDs%dUs%dReg', dsRate, usRate);
            brainDsUsRegMask = sprintf('brainDs%dUs%dRegMask', dsRate, usRate);
            brainDsUsRegMat = sprintf('brainDs%dUs%dRegMat', dsRate, dsRate); % use the original Ds5Us5 !
            brainIso2DsUssizeReg = sprintf('brainIso2Ds%dUs%dsizeReg', dsRate, usRate);

            % apply "dsXusX" registration to modality and to mask
            if usRate == usRatesSorted(1)
                md.register(brainDsUs, atlfile, 'rigid', 'multimodal', ...
                    'saveModality', brainDsUsReg, 'savetformModality', brainDsUsRegMat);
            else
                md.register(brainDsUs, atlfile, 'rigid', 'multimodal', ...
                    'saveModality', brainDsUsReg, 'loadtformModality', brainDsUsRegMat);
            end
            md.register(brainDsUsMark, atlfile, 'rigid', 'multimodal', ...
                'saveModality', brainDsUsRegMask, 'loadtformModality', brainDsUsRegMat);
            md.register(brainIsoDsUssize, atlfile, 'rigid', 'multimodal', ...
                'saveModality', brainIso2DsUssizeReg, 'loadtformModality', brainDsUsRegMat);
        end
        % view3Dopt(md.loadVolume(brainDsUs, 1), md.loadVolume(brainDsUsMask, 1));
        % view3Dopt(md.loadVolume(brainDsUsReg, 1), md.loadVolume('brainDs5Us2RegMask', 1))
    end
    
    %% seg
    if ismember('seg-inprogress', steps); % unfinished.
        for i = 1:md.getNumSubjects()
            nii = md.loadModality('seg', i);
            vol = nii.img(bbrange{i}{:});
            md.saveModality(makeNiiLike(vol, nii), 'procBrainSeg', i);
        end
        
        % crop iso to match DsXUsX
        pni = @(x, y, m, v) padNii(x, y, m, size(nii2vol(v)) - size(nii2vol(x)), 'post');
        brainCropped = sprintf('brainCropped%dSeg', dsRate);
        brainCroppedMask = sprintf('brainCropped%dMask', dsRate);
        md.applyfun(pni, {'procBrainSeg', brainCropped, brainCroppedMask, sprintf('brainDs%dUs%d', dsRate, dsRate)});
        
        % propagate registration
        usRate = dsRate;
        atlfile = eval(sprintf('atlmods.BUCKNER_ATLAS_BRAIN_PROC_DS%d_US%d', dsRate, usRate));
        brainDsUsRegMat = sprintf('brainDs%dUs%dRegMat', dsRate, dsRate); % use the original Ds5Us5!
        brainDsUsRegSeg = sprintf('brainIso2Ds%dUs%dsizeRegSeg', dsRate, dsRate); % use the original Ds5Us5!
        md.register(brainCropped, atlfile, 'rigid', 'multimodal', ...
            'saveModality', brainDsUsRegSeg, 'loadtformModality', brainDsUsRegMat, 'registeredVolumeInterp', 'nearest');
    end
    
%% mdInterpmatWarp(md, dsmod, dsusmaskmod, rmod, atlvol, regmod)
% Transform medicalDataset sparse-slice volumes according to sparse-slice interpolant matrix, as
% opposed to warping dsXusX images

    dsmod = 'brainDs5Us5Mask';
    dsusmaskmod = 'brainDs5Us5InterpMat';
    rmod = 'brainDs5Us5Regwcor';
    atlvol = BUCKNER_ATLAS_MODS.BUCKNER_ATLAS_BRAIN_PROC_DS5_US5;
    regmod = 'brainDs5Us5InterpReg';

    if ismember('mdInterpmatWarp', steps); % unfinished.

        atlnii = loadNii(atlvol);

        vi = verboseIter(1:md.getNumSubjects, 2);
        while vi.hasNext()
            i = vi.next(); 
            dsSubjNii = md.getModality(dsmod, i);
            dsusSubjmasknii = md.getModality(dsusmaskmod, i);
            interpSubjFile = md.getModality(rmod, i);
            regOutFile = md.getModality(regmod, i);

            paffine.warpvol(dsSubjNii, dsusSubjmasknii, interpSubjFile, ...
                atlnii, regOutFile);

            % visualize if return
            % q = md.loadVolume('brainDs5Us5Reg', i);
            % qn = q; q(isnan(warpedVol)) = nan;
            % qiso = md.loadVolume('brainIso2Ds5Us5sizeReg', i);
            % qiso(isnan(warpedVol)) = nan;
            % view3Dopt(warpedVol, q, qn, qiso);
        end

        vi.close();
    end


    %% copy some specific niftis to matfile
    % TODO - save all the modalities to the matfile? (simpler code) the problem is this can be huge.
    if ismember('matfile', steps);
        for usRate = usRates;
            vi = verboseIter(1:md.getNumSubjects);
            while vi.hasNext();
                i = vi.next();
                matfilename = md.getModality('matfile', i);
                vars = {};

                % original data stuff.
                eval(sprintf('brainCropped%dnii = md.loadModality(''brainCropped%d'', i);', dsRate, dsRate));
                eval(sprintf('brainCropped%d = brainCropped%dnii.img;', dsRate, dsRate));
                eval(sprintf('brainCropped%d_hdr = brainCropped%dnii.hdr;', dsRate, dsRate));
                vars = [vars, sprintf('brainCropped%d', dsRate), sprintf('brainCropped%d_hdr', dsRate)];

                % brain isotropic to the size of dsXusY
                eval(sprintf('brainIso2Ds%dUs%dsizenii = md.loadModality(''brainIso2Ds%dUs%dsize'', i);', dsRate, usRate, dsRate, usRate));
                eval(sprintf('brainIso2Ds%dUs%dsize = brainIso2Ds%dUs%dsizenii.img;', dsRate, usRate, dsRate, usRate));
                eval(sprintf('brainIso2Ds%dUs%dsize_hdr = brainIso2Ds%dUs%dsizenii.hdr;', dsRate, usRate, dsRate, usRate));
                vars = [vars, sprintf('brainIso2Ds%dUs%dsize', dsRate, usRate), sprintf('brainIso2Ds%dUs%dsize_hdr', dsRate, usRate)];

                eval(sprintf('bbrainCropped%d_sigma2 = volblur(brainCropped%d, 2);', dsRate, dsRate));
                vars = [vars, sprintf('bbrainCropped%d_sigma2', dsRate)];

                eval(sprintf('brainDs%dUs%dMasknii = md.loadModality(''brainDs%dUs%dMask'', i);', dsRate, usRate, dsRate, usRate));
                eval(sprintf('brainDs%dUs%dMask = brainDs%dUs%dMasknii.img;', dsRate, usRate, dsRate, usRate));
                eval(sprintf('brainDs%dUs%dMask_hdr = brainDs%dUs%dMasknii.hdr;', dsRate, usRate, dsRate, usRate));
                vars = [vars, sprintf('brainDs%dUs%dMask', dsRate, usRate), sprintf('brainDs%dUs%dMask_hdr', dsRate, usRate)]; 

                % registered stuff
%                 rigidRegBrainnii = md.loadModality('rigidRegBrain', i);
%                 rigidRegBrain = rigidRegBrainnii.img;
%                 rigidRegBrain_hdr = rigidRegBrainnii.hdr;
%                 vars = [vars, 'rigidRegBrain', 'rigidRegBrain_hdr'];
% 
%                 brigidRegBrain_sigma2 = volblur(rigidRegBrain, 2);
%                 vars = [vars, 'brigidRegBrain_sigma2', 'rigidRegBrain_hdr'];

                eval(sprintf('brainDs%dUs%dRegnii = md.loadModality(''brainDs%dUs%dReg'', i);', dsRate, usRate, dsRate, usRate));
                eval(sprintf('brainDs%dUs%dReg = brainDs%dUs%dRegnii.img;', dsRate, usRate, dsRate, usRate));
                eval(sprintf('brainDs%dUs%dReg_hdr = brainDs%dUs%dRegnii.hdr;', dsRate, usRate, dsRate, usRate));
                vars = [vars, sprintf('brainDs%dUs%dReg', dsRate, usRate), sprintf('brainDs%dUs%dReg_hdr', dsRate, usRate)]; 

                eval(sprintf('brainDs%dUs%dRegMasknii = md.loadModality(''brainDs%dUs%dRegMask'', i);', dsRate, usRate, dsRate, usRate));
                eval(sprintf('brainDs%dUs%dRegMask = brainDs%dUs%dRegMasknii.img;', dsRate, usRate, dsRate, usRate));
                eval(sprintf('brainDs%dUs%dRegMask_hdr = brainDs%dUs%dRegMasknii.hdr;', dsRate, usRate, dsRate, usRate));
                vars = [vars, sprintf('brainDs%dUs%dRegMask', dsRate, usRate), sprintf('brainDs%dUs%dRegMask_hdr', dsRate, usRate)];

                eval(sprintf('brainIso2Ds%dUs%dsizeRegnii = md.loadModality(''brainIso2Ds%dUs%dsizeReg'', i);', dsRate, usRate, dsRate, usRate));
                eval(sprintf('brainIso2Ds%dUs%dsizeReg = brainIso2Ds%dUs%dsizeRegnii.img;', dsRate, usRate, dsRate, usRate));
                eval(sprintf('brainIso2Ds%dUs%dsizeReg_hdr = brainIso2Ds%dUs%dsizeRegnii.hdr;', dsRate, usRate, dsRate, usRate));
                vars = [vars, sprintf('brainIso2Ds%dUs%dsizeReg', dsRate, usRate), sprintf('brainIso2Ds%dUs%dsizeReg_hdr', dsRate, usRate)];
                
                if sys.isfile(matfilename)
                    save(matfilename, vars{:}, '-append', '-v7.3');
                else
                    save(matfilename, vars{:}, '-v7.3');
                end
            end
            vi.close();
        end
    end

    %% compare results
    if ismember('visualize', steps);
        % visual test
        vol1 = md.loadVolume('brainDs5Us5', 1);
        vol2 = md.loadVolume('brainCropped5', 1);
        view3Dopt(vol1, vol2, vol2 - vol1);

        % more visualization
        atlasvol = nii2vol(BUCKNER_ATLAS_BRAIN_PROC);
        regFirstVol1 = md.loadVolume('regBrainDs5Us5', 2);
        regSecondVol1 = md.loadVolume('brainDs5Us5Reg', 2);
        regVol1 = md.loadVolume('rigidRegBrain', 2);
        regFirstMask1 = md.loadVolume('regBrainDs5Us5Mask', 2);
        regSecondMask1 = md.loadVolume('brainDs5Us5RegMask', 2);
        view3Dopt(atlasvol, regFirstVol1, regSecondVol1, regVol1, regFirstMask1, regSecondMask1);
    end
    