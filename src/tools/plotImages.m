function plotImages(data, xi, shape, varargin)
% show images on a plot
%
% data is N x D (D is dimension of images = prod(shape))
% xi is M x 2
% shape is 1 x 2
%
% Param/Value:
%   border scalar default:1
%   borderColor 1x3 vec default: [1, 0, 0] (red)
%   xscale scalar default: xspan/2.5 (xspan is maxx - minx)
%   yscale scalar default: yspan/2.5
%
% Required: addborder 
%   http://www.mathworks.com/matlabcentral/fileexchange/21005-draw-a-border-around-an-image/
    
    % parse inputs
    opts = parseInputs(data, xi, shape, varargin{:});

    % show images on the same plot.
    hold on;
    for i = 1:opts.sampleSpacing:size(xi, 1) 
        im = reshape(data(i, :), shape);
        im = imresize(im, [25, round(25./size(im, 1)*size(im, 2))], 'nearest');
        im = repmat(im, [1, 1, 3]);
        
        % add border
        if opts.border > 0
            im = addborder(im, opts.border, opts.borderColor, 'outer');
        end
        
        % get location to show image at
        xloc = xi(i, 1) + opts.xscale * [-1, 1];
        yloc = xi(i, 2) + opts.yscale * [-1, 1];
        
        % display image
        image(xloc, yloc, im); 
    end
    
    % force gray colormap
    colormap gray; 
end

function opts = parseInputs(data, xi, shape, varargin)

    p = inputParser();
    p.addRequired('data', @ismatrix);
    p.addRequired('xi', @(x) ismatrix(x) && size(x, 2) == 2);
    p.addRequired('shape', @(n) numel(n) == 2);
    p.addParameter('nSpace', [10, 10], @isvector);
    p.addParameter('border', 0, @isscalar);
    p.addParameter('sampleSpacing', 1, @isscalar);
    p.addParameter('xscale', [], @isscalar);
    p.addParameter('yscale', [], @isscalar);
    p.addParameter('borderColor', [1, 0, 0], @isvector);
    p.parse(data, xi, shape, varargin{:});
    opts = p.Results;
    
    % compute scales 
    if isempty(opts.xscale) || isempty(opts.yscale)
        % opts.nSpace = [numel(unique(xi(:, 1))), numel(unique(xi(:, 2)))];
        scales = (max(xi) - min(xi)) ./ (opts.nSpace * 2 * 1.5);
        opts.xscale = ifelse(isempty(opts.xscale), scales(1), opts.xscale);
        opts.yscale = ifelse(isempty(opts.yscale), scales(2), opts.yscale);
    end
    
end
