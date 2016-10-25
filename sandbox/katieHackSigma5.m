nel = 125; 
tmp = rand(nel, nel); 
sigma = tmp'*tmp; 
sigma = gmmIso.sigma(:, :, 1);

nTypes = 5000;
nSamp = 1; % vary this number 1, 10, 100

covHack = zeros(nel,nel); 
covTrue = zeros(nel,nel);  
covNoisy = zeros(nel,nel);  
meanD = zeros(nel, nel);
for i=1:nTypes
    i/nTypes
    %tmp = 1*rand(nel,nel); % vary this scalar to add more noise
    %D = tmp'*tmp; 
    w = rand(nel, 1);
    w(w < 0.1) = 0.1;
    D = fn3(w);
    meanD = meanD + D;
    
    covLocal = zeros(nel,nel); 
    for j=1:nSamp
        
        x = mvnrnd(zeros(nel,1), sigma); 
        n = mvnrnd(zeros(nel,1), D); 
        y = x + n; 
        
        
        covLocal = covLocal + y'*y; 
        covTrue = covTrue + x'*x; 

    end
    covLocal = covLocal/nSamp; 
    
    covHack = covHack + (covLocal - D); 
    covNoisy = covNoisy + covLocal; 
    
end
meanD = meanD ./ nTypes;

covHack = covHack./nTypes; 
covNoisy = covNoisy./nTypes; 
covTrue = covTrue./(nTypes*nSamp); 

figure;imagesc([sigma covNoisy covTrue covHack]); 
title('original | estimated from noisy y | estimated from clean x | estimated by - D hack'); 

figure;imagesc([sigma sigma+meanD covTrue covHack]); 
title('original | orig+avgD | estimated from clean x | estimated by - D hack'); 
