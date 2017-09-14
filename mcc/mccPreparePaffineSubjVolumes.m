function mccPreparePaffineSubjVolumes(varargin)
% matlab executable wrapper for preparing paffine volumes

    tic
    paffine.prepSubjVolumes(varargin{:});    
    printf('done mccPreparePaffineSubjVolumes in %3.2f\n', toc);
