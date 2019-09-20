function [ outliers, k ] = pruningNode( parent, profile )
%PRUNINGNODE: all nodes with no neighbours are deleted from the tree
%   If a node has no neighbours it means that they had only too old
%   connection that has just been pruned or that they has won only once and
%   so they coincide with a datum. In this situation the node is deleted
%   from the tree and all the points assigned are labelled as potential
%   outliers
    flag = 1;
    outliers = [];
    k = 0;
    while flag == 1
        % Iterate until there are no nodes deleted in an iteration
        flag = 0;
        % Sort children by cardinality.
        Children = parent.Children;
        childrenCard = zeros(numel(Children),1);
        for h = 1:numel(Children)
            childrenCard(h) = Children(h).VoronoiCard;
        end
        [~, idx] = sort(childrenCard, 'ascend');
        
        for i = 1:numel(Children)
            child = Children(idx(i));
            if child.NeighboursSize <= 0 || child.VoronoiCard <= profile 
    %             if child.VoronoiCard <= 1
    %                 fprintf("New outlier\n");
    %             end
%                 assert(size(parent.Children, 1) > 0)
                pOutliers = child.VoronoiSet(1 : child.VoronoiCard);
                
                for j = 1 : child.VoronoiCard
                    point = child.VoronoiSet(j);
                    point.deleteNode();
                end
                if child.NeighboursSize ~= 0 
                    for j = 1 : child.NeighboursSize
                        child.Edges(1).deleteNode();
                    end
                end
                parent.deleteChild(child);
                realOut = pruningData(pOutliers, parent, parent.Root);
                outliers = [outliers; realOut];
                k = k+1;
                flag = 1;
            end
        end
    end
    

end

