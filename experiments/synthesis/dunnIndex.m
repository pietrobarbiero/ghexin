function [ dunnIndex ] = dunnIndex( leaves )
%DUNNINDEX Summary of this function goes here
%   Detailed explanation goes here
    
    
    intraclusterDistances = zeros(numel(leaves),1);
    for i = 1: numel(leaves)
        intraclusterDistancesCurrent = zeros(leaves(i).VoronoiCard);
        for j = 1: leaves(i).VoronoiCard -1
            for k = i:leaves(i).VoronoiCard
                intraclusterDistancesCurrent(i,j) = leaves(i).VoronoiSet(j).distance(leaves(i).VoronoiSet(k));
            end
        end
        intraclusterDistances(i) = max(max(intraclusterDistancesCurrent));
    end
    intraclusterMax = max(intraclusterDistances);
    
    interclusterDistances = zeros(numel(leaves));
    for i = 1: numel(leaves)-1
        for j = i : numel(leaves)
            interclusterDistances(i,j) = leaves(i).RefVector.distance(leaves(j).RefVector);
        end
    end
    interclusterMin = min(min(interclusterDistances);
    
    dunnIndex = interclusterMin/intraclusterMax;

end

