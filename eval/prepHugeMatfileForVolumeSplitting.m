function prepHugeMatfileForVolumeSplitting(mdpath, mod, outmatfile)
% mdpath = [SYNTHESIS_DATA_PATH, '/ADNI_T1_baselines/md/adalca_wholevol_restor_md_2016_03_06.mat']
% mod = 'Ds5Us5RegMask'
% outmatile = '/data/vision/polina/scratch/adalca/tmp/ADNI_T1_baselines_wholevol_Ds5Us5RegMask_volumes.mat'

    load(mdpath);
    
    % get subjects
    tic;
    nSubjects = md.getNumSubjects();
    vi = verboseIter(1:nSubjects, 2);
    while vi.hasNext()
        i = vi.next();
        if i == 1
            volume = md.loadVolume(mod, i);
            volumes = zeros([size(volume), nSubjects]);
            volumes(:,:,:,i) = volume;
        else
            volumes(:,:,:,i) = md.loadVolume(mod, i);
        end
    end
    vi.close();
    toc;
        
    % save matfile
    tic;
    save(outmatfile, 'volumes', '-v7.3');
    toc;
    