 % assumes we have some data ready in X
 % make sure that we have similar behaviour (at least similar
 clc
 
 regval = 1e-7;
 tol = 0.001;
 reps = 10;
 K = 10;
 
 g1 = fitgmdist(X, K, 'replicates', reps, 'regularizationValue', regval, 'Options', struct('Display', 'iter', 'TolFun', tol));
   
 g2 = wgmmfit(X, X*0+1, K, 'replicates', reps, 'TolFun', tol, 'wgmmfields', ...
     struct('muUpdateMethod', 'model0', 'logpUpdateMethod', 'model0', 'covarUpdateMethod', 'model0', 'sigmareg', regval));

 
 % visualizat mus
 compareGMMmeans({g1, g2}, [1, 1, size(X, 2)], {'gmm', 'wgmm'});
 view2D({g1.mu, g2.mu});
 view2D([dimsplit(3, g1.Sigma), dimsplit(3, g2.sigma)]);