function gmm = gmdist2wgmm(gmdist)
% 
%     gmsigmainv = zeros(size(gmdist.Sigma));
%     for i = 1:size(gmdist.Sigma, 3)
%         gmsigmainv(:,:,i) = inv(gmdist.Sigma(:,:,i));
%     end


    gmm = wgmm();
    gmm.params = struct('mu', gmdist.mu, 'sigma', gmdist.Sigma, 'pi', gmdist.ComponentProportion); %, 'sigmainv', gmsigmainv);
end