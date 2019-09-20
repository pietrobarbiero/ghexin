function [ fig ] = plotTreeGHNG( Model)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

        k = 1;
        Model.GraphIndex = k;
        queue{1} = Model;
        
        label(k) = size(Model.Samples,2);
        
        
       while (~isempty(queue))
           currentModel = queue{1};
           if (length(queue) > 1)
               queue(1) = []; %pop neuron
           else
               queue = []; %pop neuron
           end
           
           NdxValidNeurons = find(isfinite(currentModel.Means(1,:)));
           fatherIndex = currentModel.GraphIndex;
           
           
        for NdxNeuro = NdxValidNeurons,        
            childIndex = k + 1;
             if ~isempty(currentModel.Child{NdxNeuro})
                 currentModel.Child{NdxNeuro}.GraphIndex = childIndex;
                 label(childIndex) = size(currentModel.Child{NdxNeuro}.Samples,2); 
                 if (~isempty(queue))
                    queue = [queue, currentModel.Child(NdxNeuro)];
                 else
                    queue = [currentModel.Child(NdxNeuro)];
                 end
             else
                 label(childIndex) = sum(currentModel.Winners == NdxNeuro);
             end
             
                s(k) = fatherIndex;
                t(k) = childIndex;
                k = k + 1;
                
        end
       
       end
  
        s = s(~isnan(s));    
        t = t(~isnan(t));
%         weights = weights(~isnan(weights));
%         weights = [sampleNumber; weights];
        G = digraph(s, t);
        fig = figure;
%         if plotWeights == 1
%             plot(G,'NodeLabel',weights, 'Layout', 'layered')
%         else
            plot(G, 'Layout', 'layered', 'NodeLabel',[label])% le label dei nodi sono le cardinalità dei Voronoi dei nodi
            set(gca,'box','off');
            set(gca,'xcolor','w','ycolor','w','xtick',[],'ytick',[]);
            axis off;
end

