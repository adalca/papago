%% test embedding spaces

%% setup
if ispc % normal people running pc
    load('C:\Users\adalca\dropbox\research\patchSynthesis\code\subspace\data\patchesIso.mat')
    addpath(genpath('C:\Users\adalca\Dropbox\MATLAB\toolboxes\'))
    addpath(genpath('C:\Users\adalca\Dropbox (Personal)\MATLAB\external_toolboxes\drtoolbox'));
else % katie.
    load('/Users/klbouman/Dropbox (MIT)/subspace/data/patchesIso.mat')
    addpath(genpath('/Users/klbouman/Desktop/Projects/medicalInpainting/'));
end

%% isomap and visualization
selpatches = patches(251:end, :);

% isomap
[mappedX, mapping] = isomap(selpatches, 2, 10);
Y = log(bsxfun(@minus, mappedX', min(mappedX', [], 2)) + 1000);

% plot the scores in 2D
figuresc; plot(Y(1, :), Y(2, :), '+');  hold on; axis equal;
plotImages(selpatches(:, 1:25), Y', [5, 5], 'sampleSpacing', 50); 
spangrid = showImageSpace(mapping.X(:, 1:25), Y', [15, 15], [5, 5], 'lambda', 10);
figuresc(); imagesc(spangrid); colormap gray;

% plot 1D for first dimention
figuresc; plot(Y(1, :), Y(1, :), '+');  hold on; axis equal;
plotImages(selpatches(:, 1:25), [Y(1, :)', Y(1, :)'], [5, 5], 'sampleSpacing', 50); 
spangrid = showImageSpace(mapping.X(:, 1:25), Y(1, :)', 15, [5, 5], 'lambda', 10);

% video view
% [~, inds] = sort(Y(1,:)); 
% for i = inds(1:50:end), pause(0.1); 
%     im = whiten(reshape(mapping.X(i, 1:25), [5 5]));
%     imagesc(im, [0 1]);
% end

%% Gaussian Mixture and Visualization moved to sandboxGMM.m


%% some random katie analyses :D
return;
patches_reshape = permute(patches, [2 1]);
patches_reshape = reshape(patches, [5 5 5 4750]);

%[coeff,score,latent] = pca(patches);
[score] = pca(patches);
figure(1);plot(score(:,1), score(:,2), '+'); hold on; axis equal; 
figure(2);plot(score(:,1), score(:,2), '+'); hold on; axis equal; 
figure(3);plot(score(:,1), score(:,2), '+'); hold on; axis equal; 

randInd = ceil(rand(1,100)*size(score,1));
for i=randInd
    im1 = patches_reshape(:,:,1,i); 
    im2 = squeeze(patches_reshape(:,3,:,i)); 
    im3 = squeeze(patches_reshape(3,:,:,i)); 
    
    figure(1); image([score(i,1)-0.05 score(i,1)+0.05],[score(i,2)-0.05 score(i,2)+0.05],im1*60); colormap gray;
    %figure(2); image([score(i,1)-0.05 score(i,1)+0.05],[score(i,2)-0.05 score(i,2)+0.05],im2*60); colormap gray;
    %figure(3); image([score(i,1)-0.05 score(i,1)+0.05],[score(i,2)-0.05 score(i,2)+0.05],im3*60); colormap gray; 
    %pause(0.1); 
end

% 
% for i=1:size(patches,1)
%     im = [patches_reshape(:,:,3,inds(i)) squeeze(patches_reshape(:,3,:,inds(i))) squeeze(patches_reshape(3,:,:,inds(i)))]; 
%     
%     imagesc(im); colormap gray; pause(0.1);
% end
% Y = lle(patches',12,2);
    

%% =====================

% clear images
% count = 1; 
% img = imresize(double(rgb2gray(imread('peppers.png'))), [40 40]);    
% for i=1:size(img,1)
%     %for j=1:size(img,2)
%         im = circshift(img, [i j]); 
%         images(count,:) = im(:);
%         count = count + 1; 
%     %end
% end
for theta = 1:270
    images(theta,1) = sin(theta*pi/180);
    images(theta,2) = cos(theta*pi/180); 
    images(theta,3) = 0; randn*0.00001; 
end


indices = randperm(99); 
[~, indicesback] = sort(indices);
origvec = 1:99;
indicesback2 = origvec(indices); 

clear images
listing = dir('/Users/klbouman/Desktop/Projects/medicalInpainting/data/teststatue/*.tiff');
for i=1:length(listing)
    im = double(imread(sprintf('/Users/klbouman/Desktop/Projects/medicalInpainting/data/teststatue/%s', listing(i).name)));
    im = imresize((im(:,:,1)), [100 100]);
    images(i,:) = im(:);
end

[score] = pca(images,2);
%Y = score'; 

%score = score(end-70:end,:);
%images = images(end-70:end,:);

images = images(indices,:);

Y=lle(images',10,1);


[mappedImages, mapping] = isomap(images, 2, 7);
Y = mappedImages';

%Y(2,:) = 1; 

[~, inds] = sort(Y(1,:)); 

clf; subplot(221);plot(Y(1,:), Y(2,:), '+'); axis equal; hold on;
rowimgs = []; 
for j=1:1:length(inds)
    i = inds(j);
    im = reshape(images(i,:), [100 100]); 
    
    subplot(222); imagesc(im); colormap gray; axis equal; axis off; title(Y(1,i)); 
    subplot(221);image([Y(1,i)-0.1 Y(1,i)+0.1],randn+[Y(2,i)-0.1 Y(2,i)+0.1],flipud(im)./255*60); axis equal;
    
    rowimgs = [rowimgs im];
    subplot(212);imagesc(rowimgs); axis equal; axis off; 
    pause(0.1); 
end

scale1 = 1e4; %1e-8; 
figure;plot(Y(1,indicesback)+1, '-*'); hold on;
for j = 1:4:numel(indicesback)
    k = indicesback(j);
    im = reshape(images(k,:), [100 100]); 
    image([j-1 j+1],([Y(1,k)-scale1*0.1 Y(1,k)+scale1*0.1])+1,flipud(im)./255*60); 
end

