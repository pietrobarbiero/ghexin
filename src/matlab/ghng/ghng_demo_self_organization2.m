% GHNG demo for self-organization with 'X' letter shape

clear all
NumSamples=312;
MaxNeurons = 20; % Maximum number of neurons in each graph
Tau = 0.1;

% The following values of the parameters are those considered in the
% original GNG paper by Fritzke (1995)
Lambda=100;
Epochs=15;
EpsilonB=0.1;
EpsilonN=0.01;
Alpha=0.5;
AMax=50;
D=0.995;

% Generate data ('X' letter shape)

load('spiral.mat')
M = X';
idx = randperm(size(M,2));    
Samples = M(:, idx);
Samples1 = GenerateSamplesImg('X_letter.bmp',NumSamples);
Samples = mapminmax(Samples,0,1);

% GHNG Training
[Model] = TrainGHNG(Samples1,Epochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);

% Plot the Model
MyAxis = [0 1 0 1];
Handle = PlotGHNG(Model,MyAxis);
h = plot(Samples1(1,:),Samples1(2,:),'*y');
uistack(h,'bottom');
axis([0 1 0 1])
hold off

