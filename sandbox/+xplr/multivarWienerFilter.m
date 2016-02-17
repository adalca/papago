%% Testing a multivariate wiener filter
% The model for our observation is:
% we have some Xt drawn from a gaussian with mean mu and sigma S, 
% and Xo is drawn from Xt with via D

% setup
nSamples = 1000; 
nDims = 100; 

%% simulate data

% first gaussian parameters
sigma = rand(nDims,nDims);
sigma = sigma*sigma'; 
mu = randn(nDims,1); 

% simulate Xt
Xt = mvnrnd(mu,sigma, nSamples); 

% simulate W.
w = rand(nDims, 1); 

% simulate observed
D = diag((-log(w)).^2); 
Xo = mvnrnd(Xt, D); 

%% test wiener filter
% Xfilt = sigma * inv(sigma + D) * (x0 - mu) + mu;
Xfilt = bsxfun(@plus, sigma / (sigma + D) * (bsxfun(@minus, Xo', mu)), mu)'; 

% see the results
plot(repmat(w', [nSamples 1]), (Xt - Xfilt).^2, '.'); 

figure();
subplot(121); plot(sum((Xt - Xo).^2,1)/nSamples, '.'); title('original error');
subplot(122); plot(sum((Xt - Xfilt).^2,1)/nSamples, '.'); title('filtered error');

%% missing data test: weights are either 1 or ~0.

% simulate weights
w = ones(nDims, 1);
idx = randsample(nDims, round(nDims)/2); 
idx = 10;
w(idx) = 0.0001;

% simulate observed
D = diag((-log(w)).^2); 
Xo = mvnrnd(Xt, D); 

% filter
Xfilt = bsxfun(@plus, sigma / (sigma + D)*(bsxfun(@minus, Xo', mu)), mu)'; 

% plot the difference with original and with 
figure();
subplot(121); plot(sum((Xt - Xo).^2,1)/nSamples, '.'); title('original error');
subplot(122); plot(sum((Xt - Xfilt).^2,1)/nSamples, '.'); title('filtered error');


%%

newsigma = sigma*inv(sigma + D)*D;
figure;plot(exp(-sqrt(diag(D))), exp(-sqrt(diag(newsigma))), '+')





