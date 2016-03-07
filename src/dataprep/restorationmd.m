function md = restorationmd(dsAmounts, buildpath, savepath, name)
% prepare a medicalDataset object for the restoration project. 
%   md = restorationmd(dsAmounts, buildpath) prepare and build a medicalDataset
%   object for the restoration project. 
%
%   md = restorationmd(dsAmounts, buildpath, savepath, name) also save the md
%   object.
%
%   md = restorationmd(2:7, [SYNTHESIS_DATA_PATH, '/buckner/proc/wholevol'], [SYNTHESIS_DATA_PATH, '/buckner'], 'wholevol');
%
%   dsAmounts - a vector of the downsampling amounts. e.g. 2:5
%   buildpath e.g.: ADNI_PATH_PROC
%   savepath e.g.: SYNTHESIS_DATA_PATH
%   name - the name of the md, mostly for saving purposes
%
% Image Restoration project

    %% create adni medicalDatasets for data.

    % initialize basic modalities
    md = medicalDataset();    
    
    % volumes
    md.addRequiredModality('orig', '%s.nii.gz');
    md.addModality('proc', '%s_proc.nii.gz');
    
    md.addModality('seg', '%s_seg.nii.gz');
    md.addModality('procSeg', '%s_proc_seg.nii.gz');
    
    % modalities matfile
    md.addModality('matfile', '%s_vols.mat');

    %% downsampled
    for s = dsAmounts % downsample amount

        % ds data in z direction
        mod = sprintf('%s_proc_ds%d.nii.gz', '%s', s); 
        md.addModality(sprintf('Ds%d', s), mod);
        
        mod = sprintf('%s_proc_ds%d_seg.nii.gz', '%s', s); 
        md.addModality(sprintf('Ds%dSeg', s), mod);

        % cropped iso to DsXUsX
        mod = sprintf('%s_cropped%d.nii.gz', '%s', s);
        md.addModality(sprintf('cropped%d', s), mod);

        % segmentation
        mod = sprintf('%s_cropped%d_seg.nii.gz', '%s', s);
        md.addModality(sprintf('cropped%dSeg', s), mod);            
        
        % cropping mask
        mod = sprintf('%s_cropped%d_cropmask.nii.gz', '%s', s);
        md.addModality(sprintf('cropped%dMask', s), mod);
    end
    
    %% ds us 
    for s = dsAmounts % downsample amount
        for u = 1:s % upsample amount
            % in most cases, all volumes will be ds by a factor of 1x1xs, then "upsampled" 
            % by a factor of 1x1xu
            %
            % thus, if s = u, you get volumes that are the same size as you started.
            %
            % exception are the *iso* volumes, which are ds by a factor of 1x1xs, then
            % ds by uxuxs then upsampled to uxuxu. 

            mod = sprintf('%s_iso_2_ds%d_us%d_size.nii.gz', '%s', s, u);
            md.addModality(sprintf('Iso2Ds%dUs%dsize', s, u), mod); 

            % this is isotropic in slices, but downsampled in z.
            mod = sprintf('%s_ds%d_isoslices%d.nii.gz', '%s', s, u); 
            md.addModality(sprintf('Ds%dIso%d', s, u), mod); 

            mod = sprintf('%s_ds%d_us%d.nii.gz', '%s', s, u);
            md.addModality(sprintf('Ds%dUs%d', s, u), mod);

            mod = sprintf('%s_ds%d_us%d_nn.nii.gz', '%s', s, u);
            md.addModality(sprintf('Ds%dUs%dNN', s, u), mod);

            mod = sprintf('%s_ds%d_us%d_dsmask.nii.gz', '%s', s, u);
            md.addModality(sprintf('Ds%dUs%dMask', s, u), mod);
            
            % segmentation
            mod = sprintf('%s_ds%d_us%d_seg.nii.gz', '%s', s, u);
            md.addModality(sprintf('Ds%dUs%dSeg', s, u), mod); 
        end
    end

    %% ds us registration modalities: ds/us and then register
    
    for s = dsAmounts % downsample amount
        for u = 1:s % upsample amount

            mod = sprintf('%s_ds%d_us%d_reg.nii.gz', '%s', s, u);
            md.addModality(sprintf('Ds%dUs%dReg', s, u), mod);

            mod = sprintf('%s_ds%d_us%d_reg.mat', '%s', s, u);
            md.addModality(sprintf('Ds%dUs%dRegMat', s, u), mod);

            mod = sprintf('%s_ds%d_us%d_dsmask_reg.nii.gz', '%s', s, u);
            md.addModality(sprintf('Ds%dUs%dRegMask', s, u), mod);

            mod = sprintf('%s_iso_2_ds%d_us%d_size_reg.nii.gz', '%s', s, u);
            md.addModality(sprintf('Iso2Ds%dUs%dsizeReg', s, u), mod); 

            mod = sprintf('%s_ds%d_us%d_regwcor.nii.gz', '%s', s, u);
            md.addModality(sprintf('Ds%dUs%dInterpReg', s, u), mod);

            mod = sprintf('%s_ds%d_us%d_reg_seg.nii.gz', '%s', s, u);
            md.addModality(sprintf('Ds%dUs%dRegSeg', s, u), mod);
        end
    end
    
    %% special ds2us2 modality 
    % Here, we downsample the iso2ds2us2size, and then ds. 
    % This is useful for testing algorithms as ds2 size, since this way the "known" 
    % planes in ds2iso match ds2us2
    
    for s = dsAmounts % downsample amount
        mod = sprintf('%s_iso_2_ds%d_us2_size__ds2us2.nii.gz', '%s', s);
        md.addModality(sprintf('Iso2Ds%dUs2size_Ds2Us2', s), mod); 

        mod = sprintf('%s_iso_2_ds%d_us2_size__ds2us2Mask.nii.gz', '%s', s);
        md.addModality(sprintf('Iso2Ds%dUs2size_Ds2Us2Mask', s), mod); 

        mod = sprintf('%s_iso_2_ds%d_us2_size__ds2us2_reg.nii.gz', '%s', s);
        md.addModality(sprintf('Iso2Ds%dUs2size_Ds2Us2Reg', s), mod); 

        mod = sprintf('%s_iso_2_ds%d_us2_size__ds2us2Mask_reg.nii.gz', '%s', s);
        md.addModality(sprintf('Iso2Ds%dUs2size_Ds2Us2MaskReg', s), mod); 
    end
    
    %% old modalities no longer in use
% 
%     % rigid registrations directly
%     md.addModality('rigidReg', '%s_reg_buckner61.nii.gz');
%     md.addModality('rigidRegBrain', '%s_brain_reg_buckner61_brain.nii.gz');
%     md.addModality('rigidRegMat', '%s_reg_buckner61.mat');
%     md.addModality('rigidRegMatBrain', '%s_brain_reg_buckner61_brain.mat');
%     
% 
%     % bounding boxes
%     md.addModality('bb', 'boundingBox/%s.nii.gz');
%     md.addModality('bbBrain', 'boundingBox/%s_brain.nii.gz');
%     md.addModality('bbmat', 'boundingBox/%s.mat');
%     md.addModality('bbmatBrain', 'boundingBox/%s_brain.mat');
%
%             % registration  modalities: first registering then propagate the ds/us
%             mod = sprintf('%s_brain_reg_ds%d_us%d.nii.gz', '%s', s, u);
%             md.addModality(sprintf('regBrainDs%dUs%d', s, u), mod);
% 
%             mod = sprintf('%s_brain_reg_ds%d_us%d_dsmask.nii.gz', '%s', s, u);
%             md.addModality(sprintf('regBrainDs%dUs%dMask', s, u), mod);
% 
%             mod = sprintf('%s_brain_reg_ds%d_us%d_nn.nii.gz', '%s', s, u);
%             md.addModality(sprintf('regBrainDs%dUs%dNN', s, u), mod);
% 
%             mod = sprintf('%s_brain_reg_ds%d_us%d_dsmask_nn.nii.gz', '%s', s, u);
%             md.addModality(sprintf('regBrainDs%dUs%dNNMask', s, u), mod);   
%
%    % add rebuilt modalities
%     for s = dsAmounts % downsample amount
%         for u = 2:s % upsample amount
%             mod = sprintf('%s_brain_ds%d_us%d_rebuilt.nii.gz', '%s', s, u);
%             md.addModality(sprintf('brainDs%dUs%d_rebuilt', s, u), mod);
% 
%             % registration modalities: ds/us and then register
%             mod = sprintf('%s_brain_ds%d_us%d_rebuilt_reg.nii.gz', '%s', s, u);
%             md.addModality(sprintf('brainDs%dUs%dReg_rebuilt', s, u), mod);
%         end
%     end
    
    %% build md
    md.build(buildpath, []);
    md.overwrite = true;
    md.verbose = 0;
    
    %% save md
    if exist('savepath', 'var')
        if ~exist('name', 'var')
            name = '';
        end
        date = datestr(now, 'yyyy_mm_dd');
        fld = [savepath, filesep, 'md', filesep];
        mkdir(fld);
        save([fld, sys.usrname, '_', name, '_restor_md_', date], 'md');
    end
    
    
    