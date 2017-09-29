function prepHugeMatfileForVolumeSplitting(mdpath, mod, outmatfile, usematfile)
% mdpath = [SYNTHESIS_DATA_PATH, '/ADNI_T1_baselines/md/adalca_wholevol_restor_md_2016_03_06.mat']
% mdpath = [SYNTHESIS_DATA_PATH, '/ADNI_T1_baselines/md/adalca_wholevol_restor_md_2016_11_12.mat']
% mod = 'Ds5Us5RegMask'
% outmatfile = '/data/vision/polina/projects/stroke/work/patchSynthesis/data/ADNI_T1_baselines/subvols/wholevol/mar12_2016/ADNI_T1_baselines_wholevol_Ds5Us5RegMask_volumes.mat'
%
% mdpath = [SYNTHESIS_DATA_PATH, '/stroke/md/adalca_brain_pad10_restor_md_2016_03_05.mat']
% mod = 'Ds7Us7Reg'
% outmatfile = '/data/vision/polina/projects/stroke/work/patchSynthesis/data/stroke/subvols/brain_pad10/mar12_2016/stroke_brain_pad10_Ds7Us7Reg_volumes.mat'

    % process inputs
    narginchk(3,4)
    if nargin < 4
        usematfile = true;
    end

    % prepare load and save files/paths
    load(mdpath);
    mkdir(fileparts(outmatfile));
    
    % get subjects
    tic;
    
    nSubjects = md.getNumSubjects();
    volIdx = find(arrayfun(@(v) sys.isfile(md.getModality(mod, v)), 1:nSubjects));
    sids = md.sids(volIdx); %#ok<NASGU>
    fprintf('found %d/%d files\n', numel(volIdx), md.getNumSubjects());
    
    % prepare large file
    volsize = size(md.loadVolume(mod, volIdx(1)));
    
    % prepare huge volumes
    if usematfile
        m = matfile(outmatfile);
        m.volumes(volsize(1), volsize(2), volsize(3), numel(volIdx)) = 0;
    else
        volumes = zeros([volsize, numel(volIdx)]);
        volumes(:,:,:,1) = volume;
    end
    
    % fill in volumes
    vi = verboseIter(volIdx, 2);
    while vi.hasNext()
        % load volumes
        [vii, idx] = vi.next();
        v = md.loadVolume(mod, vii);
        
        % add volume to huge volume
        if usematfile
            m.volumes(:, :, :, idx) = v;
        else
            volumes(:, :, :, idx) = v;
        end
    end
    vi.close();
    fprintf('Hugefile iterations took %3.2f\n', toc);
        
    % save matfile
    if ~usematfile
        tic;
        save(outmatfile, 'volumes', 'volIdx', 'sids', '-v7.3');
        fprintf('Saving hugefile took %3.2f\n', toc);
    end
    
    disp('done hugefile');
    