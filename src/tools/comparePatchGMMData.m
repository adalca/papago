function comparePatchGMMData(X1, idx1, X2, idx2, patchSize, showN)

        K = max([idx1(:);idx2(:)]);
        
        mi = min([X1(:); X2(:)]);
        ma = max([X1(:); X2(:)]);
        
        figuresc(); 
        for k = 1:K
            clstr = sprintf('mixture %d', k);

            subplot(K, 2, (k-1)*2 + 1); 
            f1 = find(idx1 == k);
            r = randperm(numel(f1));
            showN = min(showN, numel(r));
            x1 = subspacetools.reshapeN3Dto2D(X1(f1(r(1:showN)), :), patchSize);
            imagesc(x1, [mi, ma]); colormap gray; axis off equal;
            title(['samples A ', clstr]);

            subplot(K, 2, (k-1)*2 + 2); 
            f2 = find(idx2 == k);
            r = randperm(numel(f2));
            showN = min(showN, numel(r));
            x2 = subspacetools.reshapeN3Dto2D(X2(f2(r(1:showN)), :), patchSize);
            imagesc(x2, [mi, ma]); colormap gray; axis off equal;
            title(['samples B ', clstr]);
        end 
end
