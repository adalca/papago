function md = loadmd(filenames)
% filenames = [SYNTHESIS_DATA_PATH, filesep, name, filesep, sys.usrname, '_restor_md_*'];

    d = sys.fulldir(filenames);
    [~, idx] = sort(cellfun(@datenum, {d.date}), 'descend');
    latestfile = d(idx).name;

    % load md
    load(latestfile, 'md');
    