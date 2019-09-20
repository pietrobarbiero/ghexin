function [Samples]=GenerateSamplesImg(FileName,NumSamples)
% Coded by Ezequiel López-Rubio, April 2012.

Image=rot90(double(rgb2gray(imread(FileName)))/255,3);

[RowsIm,ColsIm]=size(Image);

Samples=zeros(2,NumSamples);
NumSamplesRemaining=NumSamples;
while NumSamplesRemaining>0
    PreSamples=rand(2,NumSamplesRemaining);
    Points=ceil(PreSamples.*repmat([RowsIm;ColsIm],[1 NumSamplesRemaining])); 
    LinearNdxPoints=sub2ind([RowsIm ColsIm],Points(1,:),Points(2,:));
    NdxNewSamples=find(Image(LinearNdxPoints)<1-rand(1,NumSamplesRemaining));
    NumNewSamples=numel(NdxNewSamples);
    Samples(:,NumSamples-NumSamplesRemaining+1:NumSamples-NumSamplesRemaining+NumNewSamples)=...
        PreSamples(:,NdxNewSamples);
    NumSamplesRemaining=NumSamplesRemaining-NumNewSamples;
end
        
% for Ndx=1:NumSamples
%     Hecho=0;
%     while Hecho==0
%         PreMuestra = rand(2, 1);  % Obtener NumSamples puntos de distribucion uniforme
%         Points=ceil(PreMuestra.*[RowsIm;ColsIm]);        
%         if Image(Points(1),Points(2))<(1-rand(1))
%             Samples(:,Ndx)=PreMuestra;
%             Hecho=1;
%         end
%     end
% end
