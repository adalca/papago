function gmmResultExamples(initX, compareX, resultX, trueX, gmd, errfun, patchSize, ci, dim, k)
% show top improved patch in 2D or 3D.
% ci is the cluster assignment of each patch.

    if ~exist('k', 'var'), k = 1; end
    if ~exist('dim', 'var'), dim = 2; end
    
    % compute some distances
    q1 = cellfun(errfun, dimsplit(1, compareX), dimsplit(1, trueX));
    q2 = cellfun(errfun, dimsplit(1, resultX), dimsplit(1, trueX));
    [~, mi] = sort( q1 - q2, 'descend');

    % visualize good data
    c = ci(mi(k));
    rs = @(x) reshape(x, patchSize);
    cr = @(x) x(mi(k), :);
    rscr = @(x) rs(cr(x));
    
    q = [cellfunc(rscr, {trueX, initX, compareX, resultX}), rs(gmd.mu(c,:))];
    titles = {'true', 'init', sprintf('compare (%f)', q1(k)), sprintf('result (%f)', q2(k)), 'mean'};
    if dim == 2
        
        
        figuresc();
        for i = 1:numel(q)
            subplot(5, 1, i); imagesc(subspacetools.reshape3Dto2D(q{i}(:), patchSize)); colormap gray;
            title(titles{i});
        end
        
    else
        assert(dim == 3)
        view3Dopt(q);
    end
end