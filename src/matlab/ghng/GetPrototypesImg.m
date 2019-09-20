function ImgProt = GetPrototypesImg(Centroids,Winners,ImgSize)
% Get the prototypes image from the 
% E.J. Palomo
% Inputs:
%   Centroids=Planar prototypes of the GHNG model
%   Winners=Winning neurons for each input sample
%   ImgSize=Image size
% Output:
%   ImgProt=Prototype image

% Get the prototypes of the GHNG model
Prototypes = Centroids(:,Winners);

% Prototypes redimension
ImgProt = reshape(Prototypes(1:3,:), ImgSize(3), ImgSize(1), ImgSize(2));
ImgProt = shiftdim(ImgProt,1);