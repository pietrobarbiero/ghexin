function [PSNR_ghexin, n_leaves_ghexin, ...
          PSNR_dgsot, n_leaves_dgsot, ...
          PSNR_ghng, n_leaves_ghng] = test_models_twinp()

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
rmpath(genpath('../../src/ghng'));

%% DGSOT

% clc; clear all; dbstop if error; close all;
rng(1)

addpath(genpath('../../src/dgsot'));

Samples = Samples';

classes = [];
alfa = 0.1;
sigma0 = 2;
hubHeterogenity = 0.7;
profileHeterogenity = 2;
epsilonAD = 0.05;
epsilonET = 0.03;
kappa = 1;
hetType = 1;
MaxLevels = 2;

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

MaxNeurons = 6; % Maximum number of neurons in each graph
Tau = 0.232;
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
[Model] = TrainGHNG(Samples,Epochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1,MaxLevels);

% Plot the Model
MyAxis = [0 1 0 1];
f = figure();
Handle = PlotGHNG_2(Model,MyAxis,1);
h = plot3(Samples(1,:),Samples(2,:),Samples(3,:),'*','Color',[0.9290, 0.6940, 0.1250]	);
uistack(h,'bottom');
% axis([0 1 0 1])
hold off
MyAxis = [0 1 0 1];
f = figure();
Handle = PlotGHNG_2(Model,MyAxis,2);
h = plot3(Samples(1,:),Samples(2,:),Samples(3,:),'*','Color',[0.9290, 0.6940, 0.1250]	);
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

rmpath(genpath('../../src/ghng'));
