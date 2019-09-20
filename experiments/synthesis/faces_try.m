clc; clear all; close all;

load('X_video.mat')
load('y_video.mat')

[coeff,score,latent,tsquared,explained,mu] = pca(X);

figure()
pareto(explained)
xlabel('Principal Component')
ylabel('Variance Explained (%)')

X_pca = score(:, 1:8);

save('./X_video_pca.mat', 'X_pca')