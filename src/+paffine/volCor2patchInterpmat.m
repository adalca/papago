function [R, subjInterpMask, subjPatchMins, subjPatchSize] = ...
    volCor2patchInterpmat(method, atlLoc, atlPatchSize, atlLoc2SubjSpace, varargin)
% like cor2interpmat but assumes we are given large volumes, and we want to
% do cor2interpmat on sub-patches
%
%
% TODO: sigure out how to do this without need both atl->subj *and
% subj->loc. 
%   function should just be calling src -> target ?
%   rather than worrying about bi-directionality? hard because only know
%   size a priori in atlas space


    switch method
        % obtain subj2AtlR: R subj --> atl. R is |atl|-by-|subj|
        case {'forward', 'forwardModeling', 'forwardModel'}
            % fwd
            warning('forward in currently untested!');
            [~, subjPatchMins, subjPatchSize] = paffine.atl2SubjPatch(atlLoc, atlPatchSize, atlLoc2SubjSpace);

            atlPatchRange = arrayfunc(@(a, p) a:p, atlLoc, atlLoc + atlPatchSize - 1);
            atlLoc2SubjSpacePatch = extractAndNormalizePatchCor(atlLoc2SubjSpace, atlPatchRange, subjPatchMins);
            [R, subjInterpMask, ~] = cor2interpmat(subjPatchSize, atlLoc2SubjSpacePatch);

        % obtain atl2SubjR: R atl --> subj. R is |subj|-by-|atl|
        case {'inverse', 'inverseModeling', 'inverseModel'}
            % bwd
            subjLoc2AtlSpace = varargin{1};
            [subjPatchRange, subjPatchMins, subjPatchSize] = paffine.atl2SubjPatch(atlLoc, atlPatchSize, atlLoc2SubjSpace);
            subjLoc2AtlSpacePatch = extractAndNormalizePatchCor(subjLoc2AtlSpace, subjPatchRange, atlLoc);
            [R, ~, subjInterpMask] = cor2interpmat(atlPatchSize, subjLoc2AtlSpacePatch);
     
        otherwise
            error('not a valid method for sigma estimation');
    end

end

function patchCor = extractAndNormalizePatchCor(volCor, range, mins)
    patchCor = cellfunc(@(x) x(range{:}), volCor);
    patchCor = cellfunc(@(x, m) x - m + 1, patchCor, mat2cellsplit(mins));
end
