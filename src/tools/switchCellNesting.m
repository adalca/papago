function d = switchCellNesting(c)
% switch the nesting hierarchy of a cell of cells.
% if a NxM cell of entries, each of which is 

    %
    for i = 1:numel(c)
        if ~isempty(c{i});
            for j = 1:numel(c{i})
                d{j}{i} = c{i}{j};
            end
        end
    end