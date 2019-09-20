function [ outliers ] = pruningData(data, parent, root )
%PRUNINGDATA: Re-assign potential outliers. If it fails data are outliers 
%   If a node has no neighbours it means that they had only too old
%   connection that has just been pruned or that they has won only once and
%   so they coincide with a datum. In this situation the node is deleted
%   from the tree and all the points assigned are labelled as potential
%   outliers
    outliers = [];
    for i = 1 : size(data,2)
        flag = 0;
        point = data(i);
        leaves = root.subtreeLeaves();
        
        % Try to assign data to one node of the tree 
        if size(leaves,2) == 0
           return;
        end
        dist = zeros(size(leaves, 2), 1);
        for h = 1 : size(leaves, 2)
            dist(h) = leaves(h).RefVector.distance(point);
        end
        [~,idx] = sort(dist,'ascend');
        winner1 = leaves(idx(1));
        if numel(idx) > 1
            winner2 = leaves(idx(2));
        end
        for j = 1 : size(parent.Children, 1)
           if isequal(winner1, parent.Children(j))
                % Winner is a node of the parent subnet
                if winner1.RefVector.distance(point) < winner1.T
                    winner1.add2Voronoi(point);
                    if numel(idx) > 1 && isequal(winner1.Parent, winner2.Parent)
                        winner1.addEdge(winner2);
                    end
                else
                    % The distance is greater than Threshold --> outlier
                    outliers = [outliers; point];
                    point.outlier = 1;
                end
                flag = 1;
            end
        end
        
        % Winner is a node outside parent subnet: just assign to it
        if flag == 0
            % Delete from Voronoi of the parent otherwise it will be
            % reassigned next epoch
            for k = 1 : parent.VoronoiCard
                if parent.VoronoiSet(k) == point
                   parent.deleteFromVoronoi(parent.VoronoiSet(k),k);
                   break;
                end
            end
            winner1.add2Voronoi(point);
        end
    end

end

