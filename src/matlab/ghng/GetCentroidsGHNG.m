function [NewCentroids] = GetCentroidsGHNG(Model)
% Get recursively the centroids of the GHNG model
% E.J. Palomo
% Inputs:
%   Model=GHNG model
% Output:
%   NewCentroids=Planar prototypes of the GHNG model

NdxValidNeurons = find(isfinite(Model.Means(1,:)));
NewCentroids = [];

for NdxNeuro=NdxValidNeurons,        
    if ~isempty(Model.Child{NdxNeuro}),        
        ChildCentroids = GetCentroidsGHNG(Model.Child{NdxNeuro}); 
        NewCentroids = [NewCentroids ChildCentroids];
    else
        NewCentroids = [NewCentroids Model.Means(:,NdxNeuro)];
    end    
end
