% function [stats] = ghexin_stats()
% First row: # of convex hull uses
% Second row: # of average threshold uses
% Third row: # of unexplicable points

stats = [];

%% GHEXIN

clc; clear all; dbstop if error; close all;
rng(1)

addpath(genpath('../../src/gh-exin'));
addpath(genpath('../../data'));

% Generate data ('X' letter shape)
NumSamples=312;
Samples = GenerateSamplesImg('X_letter.bmp',NumSamples)';

% GHexin parameters
HMax = 0.00002; %0.6; good for expression %0.1; good for curvatures %0.4 good for all matrix; 0.5 - not improved
Hpercentage = 0.9;
EpsilonW = 0.1;
EpsilonN = -2;
AgeMax = 5;
card = 10; %62; % 7; %10
HeightMax = 2; %1 for expression % 8 otherwise
min_epochs = 10;
avgT = true;

[nodes, leaves, ~, stats_1] = GHexinOldHcc(Samples, HMax, Hpercentage, EpsilonW, EpsilonN, AgeMax, card, HeightMax, 1, min_epochs, avgT);

nodesCard = size(nodes,2);
fig = plotGraph(nodes, nodesCard, 1);
hold off

% compute the PSNR
d = pdist(Samples, 'squaredeuclidean');
maxl2_ghexin = max(d);
% mse
mse_ghexin = 0;
Centroids_ghexin = [];
t = 1;
for i = 1 : length(leaves)
    if leaves(i).VoronoiCard > 0
        wi = leaves(i).RefVector.coordinates;
        Centroids_ghexin(:,t) = wi;
        t = t + 1;
        for j = 1 : leaves(i).VoronoiCard
            xj = leaves(i).VoronoiSet(j).coordinates;
            mse_ghexin = mse_ghexin + norm(wi - xj, 2);
        end
    end
end
mse_ghexin = mse_ghexin / NumSamples;
PSNR_ghexin = 10 * log10(maxl2_ghexin / mse_ghexin);
n_leaves_ghexin = length(Centroids_ghexin);

rmpath(genpath('../../src/gh-exin'));

%% GHEXIN

% clc; clear all; dbstop if error; close all;
rng(1)

addpath(genpath('../../src/gh-exin'));
addpath(genpath('../../data'));

% Generate data ('spiral' shape)
NumSamples=312;
load('spiral.mat')
M = X;
idx = randperm(size(M,2));
Samples = M(:, idx);
Samples = zscore(Samples);


% GHexin parameters
HMax = 0.00002; %0.6; good for expression %0.1; good for curvatures %0.4 good for all matrix; 0.5 - not improved
Hpercentage = 0.8;
EpsilonW = 0.1;
EpsilonN = -0.01; % sigma = |EpsilonN|
AgeMax = 5;
card = 10; %62; % 7; %10
HeightMax = 2; %1 for expression % 8 otherwise
min_epochs = 10;
avgT = false;

[nodes, leaves, ~, stats_2] = GHexinOldHcc(Samples, HMax, Hpercentage, EpsilonW, EpsilonN, AgeMax, card, HeightMax, 1, min_epochs, avgT);

nodesCard = size(nodes,2);
fig = plotGraph(nodes, nodesCard, 1);
hold off

% compute the PSNR
d = pdist(Samples, 'squaredeuclidean');
maxl2_ghexin = max(d);
% mse
mse_ghexin = 0;
Centroids_ghexin = [];
t = 1;
for i = 1 : length(leaves)
    if leaves(i).VoronoiCard > 0
        wi = leaves(i).RefVector.coordinates;
        Centroids_ghexin(:,t) = wi;
        t = t + 1;
        for j = 1 : leaves(i).VoronoiCard
            xj = leaves(i).VoronoiSet(j).coordinates;
            mse_ghexin = mse_ghexin + norm(wi - xj, 2);
        end
    end
end
mse_ghexin = mse_ghexin / NumSamples;
PSNR_ghexin = 10 * log10(maxl2_ghexin / mse_ghexin);
n_leaves_ghexin = length(Centroids_ghexin);

rmpath(genpath('../../src/gh-exin'));


%% GHEXIN

% clc; clear all; dbstop if error; close all;
rng(1)

addpath(genpath('../../src/gh-exin'));
addpath(genpath('../../data'));

% Generate data ('X' letter shape)
NumSamples = 1000;
Samples = betarnd(0.3, 0.3, NumSamples, 2);

% GHexin parameters
HMax = 0.00002; %0.6; good for expression %0.1; good for curvatures %0.4 good for all matrix; 0.5 - not improved
Hpercentage = 0.9;
EpsilonW = 0.35;
EpsilonN = -0.1;
AgeMax = 10;
card = 5; %62; % 7; %10
HeightMax = 2; %1 for expression % 8 otherwise
min_epochs = 10;
avgT = true;

[nodes, leaves, ~, stats_3] = GHexinOldHcc(Samples, HMax, Hpercentage, EpsilonW, EpsilonN, AgeMax, card, HeightMax, 1, min_epochs, avgT);

nodesCard = size(nodes,2);
fig = plotGraph(nodes, nodesCard, 1);
hold off

% compute the PSNR
d = pdist(Samples, 'squaredeuclidean');
maxl2_ghexin = max(d);
% mse
mse_ghexin = 0;
Centroids_ghexin = [];
t = 1;
for i = 1 : length(leaves)
    if leaves(i).VoronoiCard > 0
        wi = leaves(i).RefVector.coordinates;
        Centroids_ghexin(:,t) = wi;
        t = t + 1;
        for j = 1 : leaves(i).VoronoiCard
            xj = leaves(i).VoronoiSet(j).coordinates;
            mse_ghexin = mse_ghexin + norm(wi - xj, 2);
        end
    end
end
mse_ghexin = mse_ghexin / NumSamples;
PSNR_ghexin = 10 * log10(maxl2_ghexin / mse_ghexin);
n_leaves_ghexin = length(Centroids_ghexin);

rmpath(genpath('../../src/gh-exin'));

%% GHEXIN

% clc; clear all; dbstop if error; close all;
rng(1)

addpath(genpath('../../src/gh-exin'));
addpath(genpath('../../src/ghng'));
addpath(genpath('../../data'));

% Generate data (baboon image)
ImgOriginal = imread('baboon.tiff');
ImgOriginal = double(ImgOriginal)/255;
% ImgOriginal = ImgOriginal(128:256, 128:256, :);
Samples = reshape(shiftdim(ImgOriginal,2),3,[])';
NumSamples = size(Samples, 1);

% GHexin parameters
HMax = 0.000002; %0.6; good for expression %0.1; good for curvatures %0.4 good for all matrix; 0.5 - not improved
Hpercentage = 0.8;
EpsilonW = 0.1;
EpsilonN = -2;
AgeMax = 5;
card = 5; %62; % 7; %10
HeightMax = 2; %1 for expression % 8 otherwise
min_epochs = 10;
avgT = true;

[nodes, leaves, ~, stats_4] = GHexinOldHcc(Samples, HMax, Hpercentage, EpsilonW, EpsilonN, AgeMax, card, HeightMax, 1, min_epochs, avgT);

nodesCard = size(nodes,2);
fig = plotGraph(nodes, nodesCard, 1);
hold off

% compute the PSNR
d = pdist(Samples, 'squaredeuclidean');
maxl2_ghexin = max(d);
% mse
mse_ghexin = 0;
Centroids_ghexin = [];
t = 1;
for i = 1 : length(leaves)
    if leaves(i).VoronoiCard > 0
        wi = leaves(i).RefVector.coordinates;
        Centroids_ghexin(:,t) = wi;
        t = t + 1;
        for j = 1 : leaves(i).VoronoiCard
            xj = leaves(i).VoronoiSet(j).coordinates;
            mse_ghexin = mse_ghexin + norm(wi - xj, 2);
        end
    end
end
mse_ghexin = mse_ghexin / NumSamples;
PSNR_ghexin = 10 * log10(maxl2_ghexin / mse_ghexin);
n_leaves_ghexin = length(Centroids_ghexin);

Winners_ghexin = TestGHNG(Centroids_ghexin,Samples');
ImgProt = GetPrototypesImg(Centroids_ghexin,Winners_ghexin,size(ImgOriginal));
ImgDif = abs(ImgOriginal-ImgProt)*20;
ImgDif = imresize(ImgDif, 5);
figure, imshow(rgb2gray(ImgDif));

rmpath(genpath('../../src/gh-exin'));

%% GHEXIN

% clc; clear all; dbstop if error; close all;
rng(1)

addpath(genpath('../../src/gh-exin'));
addpath(genpath('../../src/ghng'));
addpath(genpath('../../data'));

% Generate data ('X' letter shape)
NumSamples=1000;
Samples = Generate3DSamples(5,NumSamples)';

% GHexin parameters
HMax = 0.00002; %0.6; good for expression %0.1; good for curvatures %0.4 good for all matrix; 0.5 - not improved
Hpercentage = 0.8;
EpsilonW = 0.01;
EpsilonN = 0.001;
AgeMax = 5;
card = 10; %62; % 7; %10
HeightMax = 2; %1 for expression % 8 otherwise
min_epochs = 10;
avgT = true;

[nodes, leaves, ~, stats_5] = GHexinOldHcc(Samples, HMax, Hpercentage, EpsilonW, EpsilonN, AgeMax, card, HeightMax, 1, min_epochs, avgT);

nodesCard = size(nodes,2);
fig = plotGraph(nodes, nodesCard, 1);
hold off

% compute the PSNR
d = pdist(Samples, 'squaredeuclidean');
maxl2_ghexin = max(d);
% mse
mse_ghexin = 0;
Centroids_ghexin = [];
t = 1;
for i = 1 : length(leaves)
    if leaves(i).VoronoiCard > 0
        wi = leaves(i).RefVector.coordinates;
        Centroids_ghexin(:,t) = wi;
        t = t + 1;
        for j = 1 : leaves(i).VoronoiCard
            xj = leaves(i).VoronoiSet(j).coordinates;
            mse_ghexin = mse_ghexin + norm(wi - xj, 2);
        end
    end
end
mse_ghexin = mse_ghexin / NumSamples;
PSNR_ghexin = 10 * log10(maxl2_ghexin / mse_ghexin);
n_leaves_ghexin = length(Centroids_ghexin);

rmpath(genpath('../../src/gh-exin'));
rmpath(genpath('../../src/ghng'));   

%% Stats

stats = [];
stats = [stats sum(stats_1, 2)];
stats = [stats sum(stats_2, 2)];
stats = [stats sum(stats_3, 2)];
stats = [stats sum(stats_4, 2)];
stats = [stats sum(stats_5, 2)];

c = categorical({'convex hull', 'avg-t'});
f = figure();
bar(c, stats(1:2, :))
legend('X letter','spiral','square','colors','twin-peaks')
title('Gh-exin stats - novelty')
% xlabel('experiments')
ylabel('# of events')

c = categorical({'X letter','spiral','square','colors','twin-peaks'});
f = figure();
bar(c, stats(3, :))
title('Gh-exin stats - unexplicable samples')
xlabel('experiments')
ylabel('# of events')

%% Set up

MyColor1 = [0, 0.4470, 0.7410];
MyColor2 = [0.6350, 0.0780, 0.1840];
% c = categorical({'X letter','spiral','square','colors','twin-peaks'});
% c = categorical({'convex hull', 'avg-t'});

%% X-letter stats

f = figure();
plot(stats_1(1, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor1,...
    'MarkerFaceColor',MyColor1,'MarkerSize',5,...
    'MarkerEdgeColor','none')
hold on
plot(stats_1(2, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor2,...
    'MarkerFaceColor',MyColor2,'MarkerSize',5,...
    'MarkerEdgeColor','none')
title('Gh-exin stats - X letter')
legend('convex hull', 'avg-t')
xlabel('epoch')
ylabel('# of events')
hold off

print(f, './img/ghexin_stats_Xletter_1', '-dpng');

f = figure();
plot(stats_1(3, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor1,...
    'MarkerFaceColor',MyColor1,'MarkerSize',5,...
    'MarkerEdgeColor','none')
title('Gh-exin stats - X letter')
xlabel('epoch')
ylabel('# of unexplicable samples')

print(f, './img/ghexin_stats_Xletter_2', '-dpng');

%% Spiral stats

f = figure();
plot(stats_2(1, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor1,...
    'MarkerFaceColor',MyColor1,'MarkerSize',5,...
    'MarkerEdgeColor','none')
hold on
plot(stats_2(2, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor2,...
    'MarkerFaceColor',MyColor2,'MarkerSize',5,...
    'MarkerEdgeColor','none')
title('Gh-exin stats - spiral')
legend('convex hull', 'avg-t')
xlabel('epoch')
ylabel('# of events')
hold off

print(f, './img/ghexin_stats_spiral_1', '-dpng');

f = figure();
plot(stats_2(3, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor1,...
    'MarkerFaceColor',MyColor1,'MarkerSize',5,...
    'MarkerEdgeColor','none')
title('Gh-exin stats - spiral')
xlabel('epoch')
ylabel('# of unexplicable samples')

print(f, './img/ghexin_stats_spiral_2', '-dpng');

%% Square stats

f = figure();
plot(stats_3(1, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor1,...
    'MarkerFaceColor',MyColor1,'MarkerSize',5,...
    'MarkerEdgeColor','none')
hold on
plot(stats_3(2, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor2,...
    'MarkerFaceColor',MyColor2,'MarkerSize',5,...
    'MarkerEdgeColor','none')
title('Gh-exin stats - square')
legend('convex hull', 'avg-t')
xlabel('epoch')
ylabel('# of events')
hold off

print(f, './img/ghexin_stats_square_1', '-dpng');

f = figure();
plot(stats_3(3, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor1,...
    'MarkerFaceColor',MyColor1,'MarkerSize',5,...
    'MarkerEdgeColor','none')
title('Gh-exin stats - square')
xlabel('epoch')
ylabel('# of unexplicable samples')

print(f, './img/ghexin_stats_square_2', '-dpng');

%% Colors stats

f = figure();
plot(stats_4(1, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor1,...
    'MarkerFaceColor',MyColor1,'MarkerSize',5,...
    'MarkerEdgeColor','none')
hold on
plot(stats_4(2, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor2,...
    'MarkerFaceColor',MyColor2,'MarkerSize',5,...
    'MarkerEdgeColor','none')
title('Gh-exin stats - color quantization')
legend('convex hull', 'avg-t')
xlabel('epoch')
ylabel('# of events')
hold off

print(f, './img/ghexin_stats_colors_1', '-dpng');

f = figure();
plot(stats_4(3, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor1,...
    'MarkerFaceColor',MyColor1,'MarkerSize',5,...
    'MarkerEdgeColor','none')
title('Gh-exin stats - color  quantization')
xlabel('epoch')
ylabel('# of unexplicable samples')

print(f, './img/ghexin_stats_colors_2', '-dpng');

%% Twin-peaks stats

f = figure();
plot(stats_5(1, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor1,...
    'MarkerFaceColor',MyColor1,'MarkerSize',5,...
    'MarkerEdgeColor','none')
hold on
plot(stats_5(2, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor2,...
    'MarkerFaceColor',MyColor2,'MarkerSize',5,...
    'MarkerEdgeColor','none')
title('Gh-exin stats - twin-peaks')
legend('convex hull', 'avg-t')
xlabel('epoch')
ylabel('# of events')
hold off

print(f, './img/ghexin_stats_twinp_1', '-dpng');

f = figure();
plot(stats_5(3, :),'ro-', ...
    'LineWidth',1, 'Color',MyColor1,...
    'MarkerFaceColor',MyColor1,'MarkerSize',5,...
    'MarkerEdgeColor','none')
title('Gh-exin stats - twin-peaks')
xlabel('epoch')
ylabel('# of unexplicable samples')

print(f, './img/ghexin_stats_twinp_2', '-dpng');
