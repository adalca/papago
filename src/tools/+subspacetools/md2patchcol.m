function [patches, layeridx, volidx] = md2patchcol(md, modality, patchSize, loc, pad, volnameinmatfile)
% MD2PATCHCOL extract patch column from medicalDataset object
%
% patches = md2patchcol(md, modality, patchSize, loc, pad). See detailed discussion in
% vols2patchcol(). Here, md is a medicalDataset, and modality is a string specifying which modality
% volume to work with. Note: if you have modalities that point to mat files as opposed to, say,
% niftis, this can be significantly faster. If the modality points to matfiles, the next function
% signature has to be called which specifies the volume name.
%
% patches = md2patchcol(md, modality, patchSize, loc, pad, volname) in the case that the modality
% points to a matfile, volnameinmatfile allows for the specification of the volume name in the
% matfile to be used. If the specified volume is numeric, the partial loading available through
% matfiles will allow dramatic speedup. If the specified volume is a nifti struct, unfortunately the
% partial loading does not work and the loading will be slow (as if the modality was a nifti file).
% volname can also be a cell, in which case the outputs (patches, layeridx and volidx) are each
% cells, with one for each 
%
% [patches, layeridx, volidx] = md2patchcol(...) also return layeridx and volidx (see
% vols2patchcol()).
%
% See Also: vols2patchcol, subvols2patchcol, nii2patchcol, matfile2patchcol
%
% Contact: adalca@csail

    narginchk(5, 6);

    % prepare useful numbers
    nSubjects  = md.getNumSubjects;
    nPatchesPerVol = (prod(pad * 2 + 1));
    
    % check whether modality if nifti of mat file
    [~, ~, ext] = fileparts(md.getModality(modality, 1));
    ismatfile = strcmp(ext, '.mat');
    
    % if modality is a mat file, load much faster through matfile2patchcol
    if ismatfile
        assert(nargin == 6);
        mfptrs = arrayfunc(@(x) matfile(md.getModality(modality, x)), 1:nSubjects);
        
        % check if nifti volume inside matfile
        if iscell(volnameinmatfile)
            q = whos('-file', md.getModality(modality, 1), volnameinmatfile{1});
            isnii = strcmp(q.class, 'struct');
            
            for i = 1:numel(volnameinmatfile)
                cmd = @(v) subspacetools.matfile2patchcol(mfptrs, v, patchSize, loc, pad, isnii);
                [patches{i}, layeridx{i}, volidx{i}] = cmd(volnameinmatfile);
            end
            
        else
            q = whos('-file', md.getModality(modality, 1), volnameinmatfile);
            isnii = strcmp(q.class, 'struct');
            [patches, layeridx, volidx] = ...
                subspacetools.matfile2patchcol(mfptrs, volnameinmatfile, patchSize, loc, pad, isnii);
        end
    else
        
        % prepare volumes
        patches = zeros(nSubjects * nPatchesPerVol, prod(patchSize));
        layeridx = zeros(nSubjects * nPatchesPerVol, 1);
        volidx = zeros(nSubjects * nPatchesPerVol, 1);
        
        % go over subjects, load each modality in turn, and get the patchcol. We do this in turn
        % rather than loading all the volumes to avoid memory issues.
        vi = verboseIter(1:nSubjects);
        while vi.hasNext();
            i = vi.next();
            idx = (i-1)*nPatchesPerVol+(1:nPatchesPerVol);
            
            nii = md.loadModality(modality, i);
            [patches(idx, :), layeridx(idx, :), volidx(idx, :)] = ...
                subspacetools.nii2patchcol(nii, patchSize, loc, pad);
            volidx(idx, :) = volidx(idx, :) + i - 1;
        end
        vi.close();
    end

end
