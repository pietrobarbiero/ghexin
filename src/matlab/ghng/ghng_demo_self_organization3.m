% GHNG demo for self-organization with 'Twin Peaks' 3D distribution

clear all
NumSamples=10000;
MaxNeurons = 20; % Maximum number of neurons in each graph
Tau = 0.1;

% The following values of the parameters are those considered in the
% original GNG paper by Fritzke (1995)
Lambda=100;
Epochs=2;
EpsilonB=0.2;
EpsilonN=0.006;
Alpha=0.5;
AMax=50;
D=0.995;

% Generate data ('Twin Peaks' 3D distribution)
Samples = Generate3DSamples(5,NumSamples);                      

% GHNG Training
[Model] = TrainGHNG(Samples,Epochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);

% Plot the Model
Handle = Plot3GHNG(Model);
h = plot3(Samples(1,:),Samples(2,:),Samples(3,:),'*y');               
uistack(h,'bottom');
hold off