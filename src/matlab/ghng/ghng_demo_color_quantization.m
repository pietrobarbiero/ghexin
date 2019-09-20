% GHNG demo for color quantization with the 'baboon' image

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

% Generate data ('baboon' image)
ImgOriginal = imread('baboon.tiff');
ImgOriginal = double(ImgOriginal)/255;
ImgOriginal = ImgOriginal(1:100, 1:100, :);
Samples = reshape(shiftdim(ImgOriginal,2),3,[]);

% GHNG Training
[Model] = TrainGHNG(Samples,Epochs,MaxNeurons,Tau,Lambda,EpsilonB,EpsilonN,Alpha,AMax,D,1);

% Plot the difference image amplified 20 times
Centroids = GetCentroidsGHNG(Model);
Winners = TestGHNG(Centroids,Model.Samples);
ImgProt = GetPrototypesImg(Centroids,Winners,size(ImgOriginal));
ImgDif = abs(ImgOriginal-ImgProt)*20;
figure, imshow(ImgDif);           
