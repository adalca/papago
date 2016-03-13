function vol = volresizeSimple(vol, sizeDstVol, interpMethod)
% VOLRESIZE resize N-dimentional volume
%   vol = volresize(vol, sizeDstVol) resizes the volume vol to the size sizeDstVol.
%   vol = volresize(vol, sizeDstVol, interpMethod) allows the specification of the interpolation
%       method (default is 'linear')
%       interpMethod can be anything taken by interpn(), e.g. linear or nearest, 
%       or any of these extras:
%       - 'nearestpatch' - for upsampling, use to upsample in blocks
%       - 'decimate' [TODO]: decimates every second element starting with [2, 2, 2..]
%    
%   TODO: allow for both downsampling and upsampling (in different dimensions).
%
% See Also: imBlurSep, interpn
%
% Contact: adalca@csail.mit.edu

    % check inputs
    narginchk(2, 3);
    if nargin == 2
        interpMethod = 'linear';
    end
    sizeDstVol = checkInputSizes(size(vol), sizeDstVol);

    % check that the method requires upsampling or downsampling
    % TODO - can you do both at the same time???
    sizeVol = ones(1,length(sizeDstVol));
    sizeVol(1:ndims(vol)) = size(vol); 
    isUpSampling = all(sizeVol <= sizeDstVol);
    isDownSampling = all(sizeVol >= sizeDstVol);
    assert(isUpSampling || isDownSampling, 'please only downsample or upsample the volume');
    if (isUpSampling && isDownSampling) 
        return;
    end
    isNearest = strcmp(interpMethod, 'nearest');
    isNearestPatch = strcmp(interpMethod, 'nearestpatch');
    isSimpleLinear = strcmp(interpMethod, 'simplelinear'); 
    if(isSimpleLinear)
        interpMethod = 'linear';
    end

    % if down-sampling, do a bit of smoothing first
    if isDownSampling && ~isUpSampling && ~isNearest && ~isNearestPatch && ~isSimpleLinear
        
        % blur the image. note the sigma factor is kind of arbritrary
        s = 1/4 * sizeVol ./ sizeDstVol;
        vol = volblur(vol, s);
        
        % notes from previous methods
        % s = pi^2 / 8 * sizeDstVol./sizeVol;
        % window = ceil(6 * s) + mod(ceil(6 * s), 2) + 1;
        % nn here is for the edge padding with nearest neighbours, not interpolation
    end
    
    % get the interpolation points in each dimensions
    x = cell(1, ndims(vol));
    for i = 1:ndims(vol)
        if sizeDstVol(i) > 1
            x{i} = linspace(1, size(vol,i), sizeDstVol(i));
            
            % nearest patch
            if strcmp(interpMethod, 'nearestpatch')
                assert(isUpSampling)
                ratio = sizeDstVol(i) ./ sizeVol(i);
                assert(isIntegerValue(ratio) && isodd(ratio));
                m = (ratio - 1)./2;
                
                x{i} = [ones([1, m]), linspace(1, size(vol,i), sizeDstVol(i) - 2*m), ...
                    ones([1, m]) * size(vol,i)];
            end
        else
            x{i} = (size(vol, i)+1)/2;
        end
    end
    
    % obtain a ndgrid (not meshgrid) for each dimension
    xi = cell(1, ndims(vol));
    [xi{:}] = ndgrid(x{:});

    % interpolate
    if strcmp(interpMethod, 'nearestpatch')
        interpMethod = 'nearest';
    end
    vol = interpn(vol, xi{:}, interpMethod);
    
end

function sizeDst = checkInputSizes(sizeInput, sizeDst)
    nDimsDst = numel(sizeDst);
    nDimsInput = numel(sizeInput);
    
    if nDimsDst > nDimsInput
        assert(all(sizeDst(nDimsInput + 1:end) == 1));

        msg = ['resizeNd: Destination size has more dimensions (%d) than the input vector (%d).', ...
            '\n\tSince the size in all the extra dimensions is 1, we are cropping the destination', ...
            'dimension to %d.'];
        warning('VEC:LENMATCH', msg, nDimsDst, nDimsInput, nDimsInput);
        sizeDst = sizeDst(1:nDimsInput);
    end
    
    assert(nDimsInput == numel(sizeDst));
end

%   Old Method Attempt for if (isDownSampling) code via fourier transform
%         fftVol = fftn(double(vol));
%         startVal = (floor(size(vol)/2) + 1) - ceil((sizeDstVol - 1)/2);
%         endVal = startVal + sizeDstVol - 1; 
%         fftVolCut = actionSubArray('extract', fftshift(fftVol), startVal, endVal);       
%         fftVolCutZeros = zeros(size(fftVol));
%         fftVolCutZeros = actionSubArray('insert', fftVolCutZeros, startVal, endVal, fftVolCut);
%         vol = ifftn(ifftshift(fftVolCutZeros), 'symmetric');
%         ratio = prod(sizeDstVol./size(vol));
%         vol = ifftn(ifftshift(fftVolCut), 'symmetric') .* ratio;
