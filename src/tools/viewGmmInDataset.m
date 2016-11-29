function viewGmmInDataset(atlfile, grididx, gmmfilepath, patchSize)
% 'load' or 'build'

    ax(101) = subplot(2, 4, 1); 
    ax(1) = subplot(2, 4, 1); 
    ax(2) = subplot(2, 4, 1); 
    ax(3) = subplot(2, 4, 1); 

    while true
        try
            % get the input - x, y, and mouse button used
            clear x y
            [x, y] = ginput(1);
            clickedAx = gca;
            assert(numel(x) == 1 && x > 0, 'patchview:CleanFigClose', 'unexpected input');
            assert(numel(y) == 1 && y > 0, 'patchview:CleanFigClose', 'unexpected input');
            x = round(x);
            y = round(y);
            
            % do task by axes pressed
            fidx = find(ax == clickedAx);
            switch fidx
                case 101 % control
                    % show chosen cluster. 
                    % Can also choose to show other clusters in gmm
                    % Can also show what the ISO would show?
                    % so maybe show many things? now sure.
                    showChosenCluster(loc)
                case 1
                    loc([2, 3]) = [x, y];
                case 2
                    loc([1, 3]) = [x, y];
                case 3
                    loc([1, 2]) = [x, y];
            end
            
            % draw rectangel around current pixel. 
            % TODO: should do around current patch?
            for i = 1:3
                z = loc;
                z(i) = [];
                axes(ax(i));
                patchview.drawPatchRect(z, [1, 1], 'r');
            end
             
        catch err
            okids = {'MATLAB:ginput:FigureDeletionPause', ...
                'MATLAB:ginput:FigureUnavailable', ...
                'patchview:CleanFigClose', ...
                'MATLAB:ginput:Interrupted'};
            if ~any(strcmp(err.identifier, okids))
                rethrow(err)
            end
            break;
        end
    end
end

function showChosenCluster(subvolpath, gmmpath)
    % 1. get the right gmm for this location
    
    
    
    % 2. extract gmm and find optimal logp
    % 3. dislay mu and sigma on the right. 
end
