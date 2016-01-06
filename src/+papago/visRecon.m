function visRecon(oX, oM, rX, patchSize, varargin)

    [oX, oM, rX, patchSize, nShow, dorand, titles] = parseInputs(oX, oM, rX, patchSize, varargin{:});
    
    % open new figure
    figuresc();
    
    % get order
    if dorand
        r = randperm(size(oX, 1));
    else
        r = 1:size(oX, 1);
    end
    
    % get color range
    rXc = cellfunc(@(x) x(:)', rX);
    mm = minmax([oX(:)', rXc{:}]);
    
    % get the original patches
    nX = numel(rX) + 1;
    xiso = subspacetools.reshapeN3Dto2D(oX(r(1:nShow), :), patchSize);
    subplot(2, nX, 1); 
    imagesc(xiso, mm); 
    colormap gray; axis off equal; title('original');
    
    % get the original masks
    miso = subspacetools.reshapeN3Dto2D(oM(r(1:nShow), :), patchSize);
    subplot(2, nX, nX + 1); 
    imagesc(miso, mm); 
    colormap gray; axis off equal; title('original');
    
    % go through all the rest
    for i = 1:(nX-1)
        xds = subspacetools.reshapeN3Dto2D(rX{i}(r(1:nShow), :), patchSize);
        subplot(2, nX, i+1); 
        imagesc(xds, mm); 
        colormap gray; axis off equal; 
        title(titles{i});
        
        xdel = abs(xds - xiso);
        subplot(2, nX, nX + i + 1); 
        imagesc(xdel, mm); 
        colormap gray; axis off equal; 
        title(['|orig - ', titles{i}, '|']);
    end
end

function [oX, oM, rX, patchSize, nShow, dorand, titles] = parseInputs(oX, oM, rX, patchSize, varargin)

    p = inputParser();
    p.addRequired('oX', @ismatrix);
    p.addRequired('oM', @ismatrix); % original mask
    p.addRequired('rX', @(x) ismatrix(x) || iscell(x));
    p.addRequired('patchSize', @isvector);
    p.addParameter('nShow', 10, @isscalar);
    p.addParameter('dorand', true, @islogical);
    p.addParameter('titles', {}, @(x) ischar(x) || iscellstr(x));
    p.parse(oX, oM, rX, patchSize, varargin{:});
    
    if ~iscell(rX)
        rX = {rX};
    end
    
    if isempty(p.Results.titles)
        titles = repmat({''}, [1, nX-1]);
    else
        titles = p.Results.titles;
    end
    nShow = p.Results.nShow;
    
    dorand = p.Results.dorand;
end
