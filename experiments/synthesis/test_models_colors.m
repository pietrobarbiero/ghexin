function [PSNR_ghexin, n_leaves_ghexin, ...
          PSNR_dgsot, n_leaves_dgsot, ...
          PSNR_ghng, n_leaves_ghng] = test_models_colors()

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

Winners_ghexin = TestGHNG(Centroids_ghexin,Samples');
ImgProt = GetPrototypesImg(Centroids_ghexin,Winners_ghexin,size(ImgOriginal));
ImgDif = abs(ImgOriginal-ImgProt)*20;
ImgDif = imresize(ImgDif, 5);
figure, imshow(rgb2gray(ImgDif));

rmpath(genpath('../../src/gh-exin'));   


%% DGSOT

% clc; clear all; dbstop if error; close all;
rng(1)

addpath(genpath('../../src/dgsot'));

Samples = Samples';

classes = [];
alfa = 0.2;
sigma0 = 1;
hubHeterogenity = 0.3;
profileHeterogenity = 2;
epsilonAD = 0.04;
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

Winners_dgsot = TestGHNG(Centroids_dgsot,Samples);
ImgProt = GetPrototypesImg(Centroids_dgsot,Winners_dgsot,size(ImgOriginal));
ImgDif = abs(ImgOriginal-ImgProt)*20;
ImgDif = imresize(ImgDif, 5);
figure, imshow(rgb2gray(ImgDif));

rmpath(genpath('../../src/dgsot'));

%% GHNG

% clc; clear all; dbstop if error; close all;
rng(1)

addpath(genpath('../../src/ghng'));

MaxNeurons = 7; % Maximum number of neurons in each graph
Tau = 0.1;
MaxLevels = 2;
% The following values of the parameters are those considered in the
% original GNG paper by Fritzke (1995)
Lambda=100;
Epochs=15;
EpsilonB=0.4;
EpsilonN=0.01;
Alpha=0.5;
AMax=14;
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

% Plot the difference image amplified 20 times
Centroids_ghng = GetCentroidsGHNG(Model);
Winners_ghng = TestGHNG(Centroids_ghng,Model.Samples);
ImgProt = GetPrototypesImg(Centroids_ghng,Winners_ghng,size(ImgOriginal));
ImgDif = abs(ImgOriginal-ImgProt)*20;
ImgDif = imresize(ImgDif, 5);
figure, imshow(rgb2gray(ImgDif));

rmpath(genpath('../../src/ghng'));
