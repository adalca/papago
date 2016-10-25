
X = bsxfun(@minus, bucknerIsoPatchCol, mean(bucknerIsoPatchCol, 2));
Xz = bsxfun(@minus, X, mean(X));
covar = Xz' * Xz ./ size(X, 1);
X0 = bsxfun(@minus, bucknerDsPatchCol, mean(bucknerDsPatchCol, 2));
X0z = bsxfun(@minus, X0, mean(X0));
covar0 = X0z' * X0z ./ size(X, 1);
W = bucknerDsMaskPatchCol;
covarW = W' * W ./ size(W, 1);


view2D({covar, covar0, covarW})
%%
plot(covarW(:), abs(covar(:)- covar0(:)), '.');

%% sigma3
s = 0; d = 0;
for i = 1:size(X, 1)
    x = X0z(i, :);
    ww = W(i, :)' * W(i, :);
    ww = max(ww, 0.000001);
    s = s + x(:) * x(:)' .* (ww);
    d = d + ww;
end
covar3 = s ./ d;

wtdenom = d;

%%
% note: don't use eig. it's less stable/ or doesn't work? not sure.
% TODO: exercise:
% get 100x100 covar from a 10-component PCA model
% destroy a couple of elements, and weight them low
% try to recover by doing PCA?
nComp = 41;
mic = min(covar(:));
mac = max(covar(:));

[U, S, V] = svd(covar);
recovCovar = U * S * V';

% get a smaller covar from this
ldcovar = U(:, 1:nComp) * S(1:nComp, 1:nComp) * V(:, 1:nComp)';
[U2, S2, V2] = svd(ldcovar); 
ldcovarRecon = U2(:, 1:nComp) * S2(1:nComp, 1:nComp) * V2(:, 1:nComp)';

% destory some elements
ldcovarNoisy = ldcovar;
r = rand(size(covar)); r(r>0.5) = 1;
fprintf('Destroying %d elems\n', sum(r(:)<1));
ldcovarNoisy = ldcovarNoisy .* r + (1 - r) .* 0; reshape(normrnd(ldcovarNoisy(:), std(ldcovar(:))), size(ldcovarNoisy)); 

% reconstruct via new pca
[Un, Sn, Vn] = svd(ldcovarNoisy); 
ldcovarNoisyRecon = Un(:, 1:nComp) * Sn(1:nComp, 1:nComp) * Vn(:, 1:nComp)';



ldcovarReal = covar3;
[Ur, Sr, Vr] = svd(ldcovarReal); 
ldcovarRealRecon = Ur(:, 1:nComp) * Sr(1:nComp, 1:nComp) * Vr(:, 1:nComp)';

%%
ops = {covar, recovCovar, ldcovar, ldcovarRecon, ldcovarNoisy, ldcovarNoisyRecon,  covar0, ldcovarReal, wtdenom./max(wtdenom(:)).*max(ldcovarReal(:)), ldcovarRealRecon};
ops = cellfunc(@(x) x(250:300, 1:50), ops);
titles = {'orig', 'recovOrig', 'ld', 'ldRecon', 'ldNoisy', 'ldNoisyRecon', 'covar0', 'ldcovarReal', 'wtdenom', 'ldcovarRealRecon'};
view2D(ops, 'caxis', [mic, mac], 'titles', titles); colormap gray;

%%
[U0, S0, V0] = svd(covar0);
covar0Recon = U0(:, 1:nComp) * S0(1:nComp, 1:nComp) * V0(:, 1:nComp)';

err = abs(covar-covar0); %imagesc(err./abs(covar), [0, prctile(err(:)./abs(covar(:)), 95)]);
err2 = abs(recovCovar - covar);
err3 = abs(covar0Recon - covar);

[U0, S0, V0] = svd(s./size(X, 1));
covar0Recon = U0(:, 1:nComp) * S0(1:nComp, 1:nComp) * V0(:, 1:nComp)';
err4 = abs(covar0Recon - covar);

%% from katie
ncmp = 50; 
sigma = rand(100,100);
sigma = sigma*sigma';
sigma = covar;
 
[U,S,V] = svd(sigma);
 
sigma2 = U(:,1:ncmp)*S(1:ncmp,1:ncmp)*V(:,1:ncmp)';
 
[U2, S2, V2] = svd(sigma2); 

sigma3 = U2 * S2 * V2';
 
figure; plot(diag(S(1:ncmp,1:ncmp))); hold on; plot(diag(S2(1:ncmp,1:ncmp)), 'r'); 

figure(); imagesc([sigma, sigma2, sigma3]);

%% eig via "from katie" section
ncmp = 50; 
sigma = rand(100,100);
sigma = sigma*sigma';
sigma = covar;
 
[U, S] = eig(sigma);
U = fliplr(U);
S = fliplr(flipud(S));
 
sigma2 = U(:,1:ncmp)*S(1:ncmp,1:ncmp)*U(:,1:ncmp)';
 
[U2, S2] = eig(sigma2); 
U2 = fliplr(U2);
S2 = fliplr(flipud(S2));

sigma3 = U2 * S2 * U2';

figure(); imagesc([sigma, sigma2, sigma3]);

%figure; plot(diag(S(1:ncmp,1:ncmp))); hold on; plot(diag(S2(1:ncmp,1:ncmp)), 'r'); 



