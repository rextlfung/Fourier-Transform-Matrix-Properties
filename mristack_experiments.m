%% Experiment setup

close all; clear; clc;
rng('default') % reproducibility
Nx = 256; Ny = 256; % Make small images to make things faster for testing
N_images = 21;
nbins = 100; % # of bins for histogram
%% Load in mristack example images

load mristack
images = permute(double(mristack), [2,1,3]);

% Visualize
figure(); im(images)

%% Compress each image as best rank-k approximation
k = min(Nx,Ny)/2;
images_compressed = zeros(size(images));
for n = 1:N_images
    [U,S,V] = svd(images(:,:,n));
    images_compressed(:,:,n) = U(:,1:k) * S(1:k,1:k) * V(:,1:k).';
end
images = images_compressed;

% Visualize compressed images
figure(); im(images);

%% Compute the mean, rank, Frobenius norm, 2-norm, and stable rank of each image

means = zeros(1,N_images);
ranks = zeros(1,N_images);
f_norms = zeros(1,N_images);
norms = zeros(1,N_images);
sb_ranks = zeros(1,N_images);
for n = 1:N_images
    image = double(images(:,:,n));
    means(n) = mean(image,'a');
    ranks(n) = rank(image);
    f_norms(n) = norm(image, 'fro');
    norms(n) = norm(image);
    sb_ranks(n) = (f_norms(n)/norms(n))^2;
end

% Visualize the mean pixel value of each image, and the mean rank
figure();
subplot(5,1,1); histogram(means,nbins); title(sprintf('Overall mean pixel value = %d',mean(means)));
subplot(5,1,2); histogram(f_norms,nbins); title(sprintf('Overall mean F-norm = %d',mean(f_norms)));
subplot(5,1,3); histogram(norms,nbins); title(sprintf('Overall mean 2-norm = %d',mean(norms)));
subplot(5,1,4); histogram(ranks,nbins); title(sprintf('Overall mean rank = %d',mean(ranks)));
subplot(5,1,5); histogram(sb_ranks,nbins); title(sprintf('Overall mean stable rank = %d',mean(sb_ranks)));
sgtitle('Spacial domain');
%% Compute the 2DFT of the images

kspaces = zeros(size(images));
for n = 1:N_images
    kspaces(:,:,n) = fftshift(fft2(images(:,:,n)));
end

% Visualize the k-space magnitudes
figure(); tiledlayout(ceil(sqrt(N_images)),ceil(sqrt(N_images)),'TileSpacing','none');
for n = 1:N_images
    nexttile
    im(log(abs(kspaces(:,:,n))))
    axis off;
    colormap gray;
end
%% Compute the mean, rank, Frobenius norm, 2-norm, and stable rank of each k-space

means = zeros(1,N_images);
ranks = zeros(1,N_images);
f_norms = zeros(1,N_images);
norms = zeros(1,N_images);
sb_ranks = zeros(1,N_images);
for n = 1:N_images
    kspace = double(kspaces(:,:,n));
    means(n) = mean(kspace,'a');
    ranks(n) = rank(kspace);
    f_norms(n) = norm(kspace, 'fro');
    norms(n) = norm(kspace);
    sb_ranks(n) = (f_norms(n)/norms(n))^2;
end

% Visualize the mean pixel value of each image, and the mean rank
figure();
subplot(5,1,1); histogram(abs(means),nbins); title(sprintf('Overall mean pixel magnitude = %d',mean(abs(means))));
subplot(5,1,2); histogram(f_norms,nbins); title(sprintf('Overall mean F-norm = %d',mean(f_norms)));
subplot(5,1,3); histogram(norms,nbins); title(sprintf('Overall mean 2-norm = %d',mean(norms)));
subplot(5,1,4); histogram(ranks,nbins); title(sprintf('Overall mean rank = %d',mean(ranks)));
subplot(5,1,5); histogram(sb_ranks,nbins); title(sprintf('Overall mean stable rank = %d',mean(sb_ranks)));
sgtitle('Frequency domain');