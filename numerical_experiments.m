%% Experiment setup

close all; clear; clc;
rng('default') % reproducibility
Nx = 16; Ny = 16; % Make small images to make things faster for testing
N_images = 1000;
k = 16; % Limit rank to 8 so it can increase or decrease
nbins = 200; % # of bins for histogram
%% Generate random "noise" images of specific rank

U = normrnd(0,1,[Nx, k, N_images]);
V = normrnd(0,1,[k, Ny, N_images]);
images = pagemtimes(U,V);

% Visualize the first 4 images
figure(1); montage(images(:,:,1:4)); title('First 4 images');
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
figure(2);
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

% Visualize the first 4 k-spaces
figure(3); montage(abs(kspaces(:,:,1:4))); title('Magnitude');
figure(4); montage(angle(kspaces(:,:,1:4))); title('Phase');
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
figure(5);
subplot(5,1,1); histogram(abs(means),nbins); title(sprintf('Overall mean pixel magnitude = %d',mean(abs(means))));
subplot(5,1,2); histogram(f_norms,nbins); title(sprintf('Overall mean F-norm = %d',mean(f_norms)));
subplot(5,1,3); histogram(norms,nbins); title(sprintf('Overall mean 2-norm = %d',mean(norms)));
subplot(5,1,4); histogram(ranks,nbins); title(sprintf('Overall mean rank = %d',mean(ranks)));
subplot(5,1,5); histogram(sb_ranks,nbins); title(sprintf('Overall mean stable rank = %d',mean(sb_ranks)));
sgtitle('Frequency domain');