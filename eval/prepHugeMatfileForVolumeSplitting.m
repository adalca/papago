function prepHugeMatfileForVolumeSplitting(mdpath, mod, outmatfile)
% mdpath = [SYNTHESIS_DATA_PATH, '/ADNI_T1_baselines/md/adalca_wholevol_restor_md_2016_03_06.mat']
% mod = 'Ds5Us5RegMask'
% outmatfile = '/data/vision/polina/projects/stroke/work/patchSynthesis/data/ADNI_T1_baselines/subvols/wholevol/mar12_2016/ADNI_T1_baselines_wholevol_Ds5Us5RegMask_volumes.mat'
%
% mdpath = [SYNTHESIS_DATA_PATH, '/stroke/md/adalca_brain_pad10_restor_md_2016_03_05.mat']
% mod = 'Ds7Us7Reg'
% outmatfile = '/data/vision/polina/projects/stroke/work/patchSynthesis/data/stroke/subvols/brain_pad10/mar12_2016/stroke_brain_pad10_Ds7Us7Reg_volumes.mat'

    load(mdpath);
    
    % get subjects
    tic;
    needinit = true;
    nSubjects = md.getNumSubjects();
    vi = verboseIter(1:nSubjects, 2);
    volIdx = 1:nSubjects;
    missing = false(size(volIdx));
    while vi.hasNext()
        i = vi.next();
        if ~sys.isfile(md.getModality(mod, i));
            missing(i) = true;
            continue;
        end
            
        if needinit
            volume = md.loadVolume(mod, i);
            volumes = zeros([size(volume), nSubjects]);
            volumes(:,:,:,i) = volume;
            needinit = false;
        else
            volumes(:,:,:,i) = md.loadVolume(mod, i);
        end
    end
    volIdx(missing) = []; 
    vi.close();
    toc;
    
    volumes = volumes(:,:,:,volIdx); %#ok<NASGU>
        
    % save matfile
    tic;
    mkdir(fileparts(outmatfile));
    save(outmatfile, 'volumes', 'volIdx', '-v7.3');
    toc;
    
    disp('done');
    