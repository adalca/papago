% testWgmmfit
K = 3;
D = 50;
N = 2000;
[X, mu, sigma, pi, cidx, sigmainv] = subspacetools.simgmm(N, D, K, 'simple');

% create "true wgmm"
wgmmtrue = wgmm(X, X*0 + 1, K, mu, sigma, pi, sigmainv); 

% visualize true gmm log likelihood
tll = wgmmtrue.logp();
fprintf('logp of true params: %10f\n', tll);
plot([1, 20], [tll, tll]); hold on;

%% try wgmm on simulated data with weight 1
W = ones(N, D);
[fwgmm] = wgmmfit(X, W, K, 'debug', false);
fwgmm.visualize(wgmmtrue);

% compare with matlab gmm implementation
% gmd = fitgmdist(X, k, 'Replicates', 10);

%% hide some data via weights.
onewpi = 0.2;
r = rand(N, D) < onewpi;
W = X*0;
W(r) = normrnd(1, 0.1, sum(r(:)), 1);
W(~r) = normrnd(0, 0.1, sum(~r(:)), 1);
W = within([0.1, 1], W);
% W = rand(size(X)).^3;
% W(W < 0.01) = 0;
% W = W + eps;
[fwgmm] = wgmmfit(X, W, K, 'debug', false);
fwgmm.visualize(wgmmtrue);

%% impute "missing" data
% reX = X * 0;
% mus = zeros(size(fwgmm.mu));   
% logpin = fwgmm.logpost();
% for i = 1:1;N
%     smt = 0;
%     smb = 0;
%     w = W(i, :)';
%     Dw = diag(w);
%     Dwi = diag(1./w);
%     
%     for k = 1:fwgmm.K
%         lp = logpin(i, k);
%         mm = fwgmm.sigmainv(:,:,k) * Dw * fwgmm.mu(k, :)';       
%         smt = smt + exp(lp) * mm;
%           
%         smb = smb + exp(lp) * (fwgmm.sigmainv(:,:,k));
%     end
%     reX(i, :) = Dwi * Dwi / smb * smt;
% end
% 
% X2 = W .* X + (1-W) .* reX;

warning('this doesn''t work *too* well because of incorrect sigma estimations :(')
logpin = fwgmm.logpost();
[~, mi] = max(logpin, [], 2);
reX = X * 0;
for i = 1:N
    sg = fwgmm.sigma(:,:,mi(i));
    w = W(i, :);
%     sg = diag(1./w) * sg * diag(1./w);
    xc = X(i, :) - mu(mi(i), :);
    reX(i, :) = sg * (W(i, :) .* xc)' ./ (abs(sg) * W(i, :)') + mu(mi(i), :)';
    
    % try with C' * inv(B) * b. where C are cross terms and b are all but j's element in X(i, j)
    % but somehow need to bring in weights
end
warning('bottom of division can be 0 :(');

% test

% reX = X * 0;
% for i = 1:N
% %     s = std(X(cidx == cidx(i), :), 1);
% %     sg = sigma(:,:,cidx(i)) ./ (s' * s);
%     sg = sigma(:,:,cidx(i));
%     reX(i, :) = inv(sg) * sg * X(i, :)' ;
% end 


% split
ps = dimsplit(1, p);

