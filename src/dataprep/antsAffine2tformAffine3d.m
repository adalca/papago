function tform = antsAffine2tformAffine3d(tformAntsFilename)
% Transform a 3d ANTs (ITK) affine matrix to a matlab tform T

    % read tform
    fid = fopen(tformAntsFilename);
    tline = fgetl(fid);
    T = nan;
    C = nan;
    while ischar(tline)
        disp(tline)
        str = 'Parameters:';
        if numel(tline) > numel(str) && strcmp(tline(1:numel(str)), str) == 1
            T = reshape(str2num(tline(numel(str)+1:end)), [3, 4])';
            T = [T, zeros(4, 1)];
            T(4,4) = 1;
        end
        
        str = 'FixedParameters:';
        if numel(tline) > numel(str) && strcmp(tline(1:numel(str)), str) == 1
            C = str2num(tline(numel(str)+1:end));
        end
        
        tline = fgetl(fid);
    end
    fclose(fid);
    
    
    % get offset
    newC = T(1:3,1:3) * (zeros(3, 1) - C') + C' + T(4, 1:3)';
    
    % get new transform and switch x,y as appropriate.
    newT = T;
    newT(3, 1:2) = -T([2, 1], 3);
    newT(1:2, 3) = -T(3, [2, 1]);
    newT(4, 1:3) = newC([2, 1, 3]);
    newT(4, 1:2) = -newT(4, 1:2);

    tform = affine3d(newT);