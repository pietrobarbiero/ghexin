function [ fig ] = plotGraph( nodes, sampleNumber, pointClasses )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if nargin == 2

        nodesCard = size(nodes,2);

        s = NaN(nodesCard,1);
        t = NaN(nodesCard,1);
        weights = NaN(nodesCard,1);

        k = 1;
        for i = 1 : nodesCard
            if nodes(i).IsLeaf ~= true
                for j = 1 : size(nodes(i).Children,1)
                    s(k) = nodes(i).Number;
                    t(k) = nodes(i).Children(j).Number;
%                     weights(k) = nodes(i).Children(j).VoronoiCard;
                    weights(k) = nodes(i).Children(j).Number;
                    k = k + 1;
                end
            end
        end

        s = s(~isnan(s));    
        t = t(~isnan(t));
        weights = weights(~isnan(weights));
%         weights = [sampleNumber; weights];
        weights = [1; weights];
        G = digraph(s, t);
        fig = figure;
        plot(G,'NodeLabel',weights, 'Layout', 'layered')
   
    elseif nargin == 3
        pointClasses = pointClasses + 2;
        nodesCard = size(nodes,2);

        s = NaN(nodesCard,1);
        t = NaN(nodesCard,1);
        weights = NaN(nodesCard,1);
        classes = zeros(nodesCard-1,3);
        
        k = 1;
        for i = 1 : nodesCard
            if nodes(i).IsLeaf ~= true
                for j = 1 : size(nodes(i).Children,1)
                    s(k) = nodes(i).Number;
                    t(k) = nodes(i).Children(j).Number;
                    weights(k) = nodes(i).Children(j).VoronoiCard;
                    node = nodes(i).Children(j);
                    for h = 1 : node.VoronoiCard
                        classes(k,pointClasses(node.VoronoiSet(h).number)) = classes(k,pointClasses(node.VoronoiSet(h).number)) + 1;
                    end
                    k = k + 1;
                end
            end
        end
        classes = [0,0,0; classes];
        for h = 1 : nodes(1).VoronoiCard
            classes(1,pointClasses(nodes(1).VoronoiSet(h).number)) = classes(1,pointClasses(nodes(1).VoronoiSet(h).number)) + 1;
        end

        s = s(~isnan(s));    
        t = t(~isnan(t));
        weights = weights(~isnan(weights));
        weights = [sampleNumber; weights];
        G = digraph(s, t);
        fig = figure;
        
        nodesColor = NaN(nodesCard,3);
        for i = 1 : nodesCard
            numberOfPoints = sum(classes(i,:));
            if numberOfPoints ~= 0
                nodesColor(i,:) = classes(i,:)./numberOfPoints;
            else
                nodesColor(i,:) = 0;
            end
        end
        plot(G,'NodeLabel', weights, 'Layout', 'layered', 'NodeColor', nodesColor)
    end
end

