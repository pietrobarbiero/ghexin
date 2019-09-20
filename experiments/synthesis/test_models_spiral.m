function [PSNR_ghexin, n_leaves_ghexin, ...
          PSNR_dgsot, n_leaves_dgsot, ...
          PSNR_ghng, n_leaves_ghng] = test_models_spiral()
      
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

[nodes, leaves, ~] = GHexinOldHcc(Samples, HMax, Hpercentage, EpsilonW, EpsilonN, AgeMax, card, HeightMax, 1, min_epochs, avgT);

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

%% DGSOT

% clc; clear all; dbstop if error; close all;
rng(1)
format long

addpath(genpath('../../src/dgsot'));
% addpath(genpath('../../data'));

% Generate data ('spiral' shape)
% NumSamples=312;
% load('spiral.mat')
% M = X;
% idx = randperm(size(M,2));
% Samples = M(:, idx);
Samples = Samples';

classes = [];
alfa = 0.1;
sigma0 = 0.2;
hubHeterogenity = 0.3;
profileHeterogenity = 2;
epsilonAD = 0.04;
epsilonET = 0.03;
kappa = 0;
hetType = 1;
MaxLevels = 3;

[nodes, leaves, outliers] = DGSOT(Samples, alfa, sigma0, hubHeterogenity, profileHeterogenity, epsilonAD, epsilonET, kappa, hetType, classes, MaxLevels); 

% plotGraph(nodes, nodes(1).VoronoiCard);
% hold off

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

rmpath(genpath('../../src/dgsot'));

%% GHNG

% clc; clear all; dbstop if error; close all;
rng(1)

addpath(genpath('../../src/ghng'));
% addpath(genpath('../../data'));
% 
% NumSamples=312;
% load('spiral.mat')
% M = X;
% idx = randperm(size(M,2));
% Samples = M(:, idx);
% Samples = Samples';

MaxNeurons = 35; % Maximum number of neurons in each graph
Tau = 0.001;
MaxLevels = 2;
% The following values of the parameters are those considered in the
% original GNG paper by Fritzke (1995)
Lambda=100;
Epochs=15;
EpsilonB=0.1;
EpsilonN=0.01;
Alpha=0.5;
AMax=5;
D=0.995;
% GHNG Training
[Model] = TrainGHNG(Samples,Epochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1,MaxLevels);

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

% rmpath(genpath('../../src/ghng'));
