function [ fig ] = plotGraph( nodes, sampleNumber, Reorder, pointClasses )
%PLOTGRAPH: Customized graph plot
%   It may works in 2 ways:
%% Plot graph without considering tissue classes
    if nargin == 3
        nodesCard = size(nodes,2);

        s = NaN(nodesCard,1);
        t = NaN(nodesCard,1);
        weights = strings(nodesCard);

        k = 1;
        sHighlights = [];
        tHighlights = [];
        for i = 1 : nodesCard
            if nodes(i).IsLeaf ~= true
                for j = 1 : size(nodes(i).Children,1)
                    s(k) = nodes(i).Number;
                    if Reorder == 1
                        nodes(i).Children(j).Number = k + 1;
                    end
                    if nodes(i).Highlight == 1
                        sHighlights = [sHighlights; nodes(i).Number];
                        tHighlights = [tHighlights; nodes(i).Children(j).Number];
                    end
                    t(k) = nodes(i).Children(j).Number;
%                     weights(k) = nodes(i).Children(j).VoronoiCard;
                    weights(k) = nodes(i).Children(j).Number;
%                    weights(k) = "V: "+nodes(i).Children(j).VoronoiCard+ "H: " +num2str(nodes(i).Children(j).Hcc(),3);
%                    weights(k) = "H: " + num2str(nodes(i).Children(j).Hcc(),3);
                    k = k + 1;
                end
            end
        end

        s = s(~isnan(s));    
        t = t(~isnan(t));
        weights = weights(1 : size(s,1));
        str = 1;
%         str ="V: "+nodes(1).VoronoiCard + "H: " +num2str(nodes(1).Hcc(),3);
%         str ="H: " +num2str(nodes(1).Hcc(),3);
        weights = [str weights];
        weights = cellstr(weights);
        G = digraph(s, t);
        fig = figure;
        if Reorder == 1
            p = plot(G,'NodeLabel', weights, 'Layout', 'layered');
        else
            p = plot(G, 'Layout', 'layered');
        end
        highlight(p, sHighlights, tHighlights, 'NodeColor', 'red', 'EdgeColor', 'red');
        
    %% plot nodes with colors representing the classes elseif nargin == 4
    else
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

