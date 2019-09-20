function [PSNR_ghexin, n_leaves_ghexin, ghexin_avg_silhouette, ghexin_dbindex, ET_ghexin, ...
          PSNR_dgsot, n_leaves_dgsot, dgsot_avg_silhouette, dgsot_dbindex, ET_dgsot, ...
          PSNR_ghng, n_leaves_ghng, ghng_avg_silhouette, ghng_dbindex, ET_ghng] = test_models_X(seed)

%% GHEXIN

% clc; clear all; dbstop if error; close all;
rng(seed)

addpath(genpath('../../src/matlab/gh-exin'));
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
HeightMax = 3; %1 for expression % 8 otherwise
min_epochs = 10;
avgT = true;
toPlot = false;

tic
[nodes, leaves, outliers, points, stats] = GHexinOldHcc(Samples, HMax, Hpercentage, EpsilonW, EpsilonN, AgeMax, card, HeightMax, toPlot, min_epochs, avgT);
ET_ghexin = toc

nodesCard = size(nodes,2);
fig = plotGraph(nodes, nodesCard, 1);
hold off

k = 1;
for i = 1:numel(nodes)
    if nodes(i).Height == 2 
        nodesFirstLevel(k) = nodes(i);
        k = k+1;
    end
end
plotCluster(Samples, nodesFirstLevel, 0);
hold off
plotCluster(Samples, leaves, 0);
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

% extract point-cluster membership (take into account outliers)
point_membership = get_ghexin_point_membership(points);
notNullIdx = find(not(isnan(point_membership)));
point_membership_notnull = point_membership(notNullIdx);
samples_notnull = Samples(notNullIdx,:);

% order cluster centroids according to the number
numbers = [];
for i = 1:numel(leaves)
    numbers = [numbers; leaves(i).Number];
end
[~, idx] = sort(numbers);
centroids = Centroids_ghexin(:,idx)';

% compute silhouette
figure
[ghexin_silhouette, f] = silhouette(samples_notnull, point_membership_notnull, 'Euclidean');
ghexin_avg_silhouette = mean(ghexin_silhouette);

% compute davies bouldin index
ghexin_dbindex = db_index(samples_notnull, point_membership_notnull, centroids);

rmpath(genpath('../../src/matlab/gh-exin'));

%% DGSOT

% clc; clear all; dbstop if error; close all;
rng(seed)

addpath(genpath('../../src/matlab/dgsot'));

Samples = Samples';

classes = [];
alfa = 0.2;
sigma0 = 1;
hubHeterogenity = 0.3;
profileHeterogenity = 10;
epsilonAD = 0.046;
epsilonET = 0.03;
kappa = 1;
hetType = 1;
MaxLevels = 2;

tic
[nodes, leaves, outliers, points] = DGSOT(Samples, alfa, sigma0, hubHeterogenity, profileHeterogenity, epsilonAD, epsilonET, kappa, hetType, classes, MaxLevels); 
ET_dgsot = toc

plotGraph(nodes, nodes(1).VoronoiCard);
hold off
k = 1;
clear nodesFirstLevel
for i = 1:numel(nodes)
    if nodes(i).Height == 2 
        nodesFirstLevel(k) = nodes(i);
        k = k+1;
    end
end
plotCluster(Samples, nodesFirstLevel, 0);
hold off
plotCluster(Samples, leaves, 0);
hold off

% compute the PSNR
d = pdist(Samples', 'squaredeuclidean');
maxl2_dgsot = max(d);
% mse
mse_dgsot = 0;
Centroids_dgsot = [];
t = 1;
for i = 1 : length(leaves)
    if leaves(i).VoronoiCard > 0
        wi = leaves(i).RefVector.coordinates;
        Centroids_dgsot(:,t) = wi;
        t = t + 1;
        for j = 1 : leaves(i).VoronoiCard
            xj = leaves(i).VoronoiSet(j).coordinates;
            mse_dgsot = mse_dgsot + norm(wi - xj, 2);
        end
    end
end
mse_dgsot = mse_dgsot / NumSamples;
PSNR_dgsot = 10 * log10(maxl2_dgsot / mse_dgsot);
n_leaves_dgsot = length(Centroids_dgsot);

% order cluster centroids according to the number
numbers = [];
for i = 1:numel(leaves)
    if leaves(i).VoronoiCard > 0
        numbers = [numbers; leaves(i).Number];
    end
end
[~, idx] = sort(numbers);
centroids = Centroids_dgsot(:,idx)';

% compute silhouette
point_membership = get_dgsot_point_membership(points);
figure
dgsot_silhouette = silhouette(Samples', point_membership, 'Euclidean');
dgsot_avg_silhouette = mean(dgsot_silhouette);

% compute davies bouldin index
dgsot_dbindex = db_index(Samples', point_membership, centroids);

rmpath(genpath('../../src/matlab/dgsot'));

%% GHNG

% clc; clear all; dbstop if error; close all;
rng(seed)

addpath(genpath('../../src/matlab/ghng'));

MaxNeurons = 6; % Maximum number of neurons in each graph
Tau = 0.25;
MaxLevels = 2;
% The following values of the parameters are those considered in the
% original GNG paper by Fritzke (1995)
Lambda=100;
Epochs=15;
EpsilonB=0.1;
EpsilonN=0.01;
Alpha=0.5;
AMax=50;
D=0.995;

% GHNG Training
tic
[Model] = TrainGHNG(Samples,Epochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1,MaxLevels);
ET_ghng = toc

% Plot the Model
MyAxis = [0 1 0 1];
f = figure();
Handle = PlotGHNG_2(Model,MyAxis,1);
h = plot(Samples(1,:),Samples(2,:),'*','Color',[0.9290, 0.6940, 0.1250]	);
uistack(h,'bottom');
% axis([0 1 0 1])
hold off
MyAxis = [0 1 0 1];
f = figure();
Handle = PlotGHNG_2(Model,MyAxis,2);
h = plot(Samples(1,:),Samples(2,:),'*','Color',[0.9290, 0.6940, 0.1250]	);
uistack(h,'bottom');
% axis([0 1 0 1])
hold off

% compute the PSNR
d = pdist(Samples', 'squaredeuclidean');
maxl2_ghng = max(d);
mse_ghng = 0;
[NewCentroids, mse_ghng] = GetMseGHNG(Model, mse_ghng);
mse_ghng = mse_ghng / NumSamples;
PSNR_ghng = 10 * log10(maxl2_ghng / mse_ghng);

% count the number of leaf nodes
n_leaves_ghng = length(NewCentroids);

% compute silhouette
winners = TestGHNG(NewCentroids, Samples);
figure
[ghng_silhouette, f] = silhouette(Samples', winners, 'Euclidean');
ghng_avg_silhouette = mean(ghng_silhouette);

% compute davies bouldin index
ghng_dbindex = db_index(Samples', winners, NewCentroids');

rmpath(genpath('../../src/ghng'));
