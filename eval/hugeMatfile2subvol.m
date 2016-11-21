function [subvolume, nfo] = hugeMatfile2subvol(matfilefile, subvolLoc, subvolSize, subvolfile)
% given a matfile which has the 4D variable 'volumes' (volsize x nSubjects)
% extract the subvolumes at location subvolLoc, of size subvolSize, and save the subvolumes in
% subvolfile. the subvolumes will be called 'subvolume'

    narginchk(4, 6)

    m = matfile(matfilefile);
    
    if ischar(subvolLoc), subvolLoc = str2num(subvolLoc); end
    if ischar(subvolSize), subvolSize = str2num(subvolSize); end
    
    range = arrayfunc(@(l,s) l:(l+s-1), subvolLoc, subvolSize);
    
    subvolume = m.volumes(range{:}, :);
    nfo.subvolLoc = subvolLoc;
    nfo.subvolSize = subvolSize;
    
    save(subvolfile, 'subvolume', 'nfo');
    