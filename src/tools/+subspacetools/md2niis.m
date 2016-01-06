function niis = md2niis(md, dataset, varargin)
% subspacetools.md2niis
% subspace project specific loading of niftis
%
% e.g. 
%   niis = subspacetools.md2niis(md, 'buckner', 5, 5)

    switch dataset
        case 'buckner'
            dsRate = varargin{1};
            usRate = varargin{2};

            niis.iso = cell(md.getNumSubjects, 1);
            niis.ds = cell(md.getNumSubjects, 1);
%             niis.dsnn = cell(md.getNumSubjects, 1);
            niis.mask = cell(md.getNumSubjects, 1);

            % load all the volumes
            vi = verboseIter(1:md.getNumSubjects);
            while vi.hasNext()
                [~, i] = vi.next();
                % Extract iso volume 
                niis.iso{i} = md.loadModality('rigidRegBrain', i);

                % Extract ds volume
                brainDsUsReg = sprintf('brainDs%dUs%dReg', dsRate, usRate);
                %brainDsUsReg = sprintf('regBrainDs%dUs%d', dsRate, usRate);
                niis.ds{i} = md.loadModality(brainDsUsReg, i);

                % get mask
                brainDsUsRegMask = sprintf('brainDs%dUs%dRegMask', dsRate, usRate);
                %brainDsUsRegMask = sprintf('regBrainDs%dUs%dMask', dsRate, usRate);
                niis.mask{i} = md.loadModality(brainDsUsRegMask, i);
            end

            vi.close();
            
        case 'stroke'
            error('md2niis: stroke dataset unimplemented');
            
        otherwise
            error('unknown dataset');
    end
    