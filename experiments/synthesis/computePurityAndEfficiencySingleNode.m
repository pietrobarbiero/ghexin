function [ efficiency, purity, Eidx, Pidx ] = computePurityAndEfficiencySingleNode(dataClassesLabels_Vec, voronoiSetIndices )
   %%calcola la purezza e l'efficienza di quatizzazione di un nodo della
   %%rete dati il suo Voronoi set ed il vettore delle label associate alle
   %%classi
   %param dataClassesLabels vettore grande quanto il training set.
   %this vector contains the class label of each element of the training set
   %param voronoiSetIndices sono gli indici (riga o colonna) all'interno del training set
   %dei dati del Voronoi set del nodo in esame
    
   %% init
    uniqueClasses = unique(dataClassesLabels_Vec(~isnan(dataClassesLabels_Vec)));
    numClasses = length(uniqueClasses); %sarebbe 434/7 = 62
    class = zeros(1,numClasses);
    nodeVoronoiCard = length(voronoiSetIndices);
    
    %% calcolo purezza
    j = 0;
    for h = 1 : nodeVoronoiCard
	% per ogni classe, vedo quanti punti sono presenti all'interno di questo VoronoiSet
       idx = find(uniqueClasses == dataClassesLabels_Vec(voronoiSetIndices(h)));
%        idx = dataClassesLabels_Vec(voronoiSetIndices(h));
       if ~isnan(idx)
           class(idx) = class(idx) + 1;
       else 
           j = j+1;
       end
    end
    purities = class./(nodeVoronoiCard-j); 
    
    %% calcolo efficienza
    totvalueclasses = zeros(1,numClasses);
    efficiencies = zeros(1,numClasses);

    for i = 1 : numClasses
       %how many elements of i-th class are in this neuron Voronoi set w.r.t all the elements of i-th class 
	   totvalueclasses(i) = sum(dataClassesLabels_Vec(:) == uniqueClasses(i));% sarebbe 7
        efficiencies(i) = class(i)./totvalueclasses(i);
    end
    
    %% output 
    % Calculating purity and efficiencies as the maximum purities and
    % efficiencies between classes 
    [purity, Pidx] = max(purities);
    [efficiency, Eidx] = max(efficiencies);
end

