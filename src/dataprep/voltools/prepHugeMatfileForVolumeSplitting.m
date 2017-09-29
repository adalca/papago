function prepHugeMatfileForVolumeSplitting(mdpath, mod, outmatfile)
% mdpath = [SYNTHESIS_DATA_PATH, '/ADNI_T1_baselines/md/adalca_wholevol_restor_md_2016_03_06.mat']
% mdpath = [SYNTHESIS_DATA_PATH, '/ADNI_T1_baselines/md/adalca_wholevol_restor_md_2016_11_12.mat']
% mod = 'Ds5Us5RegMask'
% outmatfile = '/data/vision/polina/projects/stroke/work/patchSynthesis/data/ADNI_T1_baselines/subvols/wholevol/mar12_2016/ADNI_T1_baselines_wholevol_Ds5Us5RegMask_volumes.mat'
%
% mdpath = [SYNTHESIS_DATA_PATH, '/stroke/md/adalca_brain_pad10_restor_md_2016_03_05.mat']
% mod = 'Ds7Us7Reg'
% outmatfile = '/data/vision/polina/projects/stroke/work/patchSynthesis/data/stroke/subvols/brain_pad10/mar12_2016/stroke_brain_pad10_Ds7Us7Reg_volumes.mat'

    load(mdpath);
    
    % get subjects
    tic;
    
    nSubjects = md.getNumSubjects();
    volIdx = find(arrayfun(@(v) sys.isfile(md.getModality(mod, v)), 1:nSubjects));
    sids = md.sids(volIdx); %#ok<NASGU>
    fprintf('found %d/%d files\n', numel(volIdx),md.getNumSubjects());
    
    volume = md.loadVolume(mod, volIdx(1));
    volumes = zeros([size(volume), numel(volIdx)]);
    volumes(:,:,:,1) = volume;
    
    vi = verboseIter(volIdx, 2);
    while vi.hasNext()
        [vii, idx] = vi.next();
        v = md.loadVolume(mod, vii);
        volumes(:,:,:,idx) = v;
    end
    vi.close();
    toc;
        
    % save matfile
    disp('saving matfile');
    
    tic;
    mkdir(fileparts(outmatfile));
    save(outmatfile, 'volumes', 'volIdx', 'sids', '-v7.3');
    toc;
    
    disp('done');
    