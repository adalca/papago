function [recovpatches, nanpatches] = sandboxDeleteReconstruct(patches, eigenvec, sigmas, nDestroy, nReps, truescores)
% given true pca (via eigenvec and sigmas), delete features and inpaint
%
% repeat nRep times:
% - delete features from patches
% - 


    recovpatches = cell(nReps, 1);
    recovscores = cell(nReps, 1);

    % repeat the experiment several times
    for r = 1:nReps

        % destroy patches
        nanpatches = destroyFeatures(patches, nDestroy, 'rand-consistent');
        invalid = isnan(nanpatches);

        % inpaint data a via p(a|b), where a is vector of missing data and b is vector of known data
        recovpatches{r} = pcax.inpaint(nanpatches', eigenvec, sigmas)';
        recoverr = bsxfun(@minus, recovpatches{r}(:, invalid{r}), patches(:, invalid{r}));

        % TODO: alternatively, try to inpaint with the means, and re-project to lower number of
        %   eigenvecs.
        
        % recover the scores as well to check
        recovscores{r} = pcax.project(recovpatches{r}', eigenvec)';
        scoreerr = bsxfun(@minus, recovscores{r}(:, 1:size(truescores, 2)), truescores);


        % display results
        figure();
        subplot(231); imagesc(invalid{r}); 
        subplot(232); imagesc(patches); colormap gray; 
        subplot(233); imagesc(recovpatches{r}); colormap gray;

        subplot(234); boxplot(scoreerr); title('Scores error');
        subplot(235); boxplot(recoverr); title('Pixel error');

        fprintf('Mean score abserror: %e\n', mean(abs(scoreerr(:))));
        fprintf('Mean pixel abserror: %e\n', mean(abs(recoverr(:))));
    end