function [NewCentroids, mse_ghng] = GetMseGHNG(Model, mse_ghng)
% Get recursively the centroids of the GHNG model
% E.J. Palomo
% Inputs:
%   Model=GHNG model
% Output:
%   NewCentroids=Planar prototypes of the GHNG model

NdxValidNeurons = find(isfinite(Model.Means(1,:)));
NewCentroids = [];

for NdxNeuro=NdxValidNeurons
    if ~isempty(Model.Child{NdxNeuro})       
        [ChildCentroids, mse_ghng] = GetMseGHNG(Model.Child{NdxNeuro}, mse_ghng); 
        NewCentroids = [NewCentroids ChildCentroids];
    else
        NewCentroids = [NewCentroids Model.Means(:,NdxNeuro)];
        wi = Model.Means(:,NdxNeuro);
        for j = 1 : length(Model.Samples)
            if Model.Winners(j) == NdxNeuro
                xj = Model.Samples(:, j);
                mse_ghng = mse_ghng + norm(wi - xj, 2);
            end
        end
    end    
end