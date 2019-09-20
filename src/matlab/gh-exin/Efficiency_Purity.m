function [ efficiency, purity, Eidx, Pidx ] = Efficiency_Purity( node, dataClasses)
   
    numClasses = size(unique(dataClasses(~isnan(dataClasses))),1);
    
    class = zeros(1,numClasses);
    j = 0;
    for h = 1 : node.VoronoiCard
       idx = dataClasses(node.VoronoiSet(h).number);
       if ~isnan(idx)
           class(idx) = class(idx) + 1;
       else 
           j = j+1;
       end
   end
    
    purities = class./(node.VoronoiCard-j); 
    totvalueclasses = zeros(1,numClasses);
    efficiencies = zeros(1,numClasses);

    for i = 1 : numClasses
       totvalueclasses(i) = sum(dataClasses(:) == i);
        efficiencies(i) = class(i)./totvalueclasses(i);
    end
    
    % Calculating purity and efficienciues as the maximum purities and
    % efficiencies between classes 
    [purity, Pidx] = max(purities);
    [efficiency, Eidx] = max(efficiencies);
  
end

