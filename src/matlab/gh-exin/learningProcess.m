function [ counter ] = learningProcess(parent, EpsilonW, EpsilonN, epoch )
%     fprintf('Learning Process\n');
    counter = 1;
    flag = 0;

    for i = 1 : parent.VoronoiCard
        if parent.VoronoiCard == 0
            fprintf('Error\n');
        end
        point = parent.VoronoiSet(i);
      
        %plotcluster - TO BE DELETED
%         if rem (i,1000) == 0
%             plotCluster(parent.voronoiMatrix, parent.Children');
%             hold on
%             coord=point.coordinates;
%             scatter(coord(1),coord(2), 50,'k','filled')
%         end
        if point.outlier == 0

            %point.coordinates
            if epoch ~= 0
                point.deleteNode();
            end
            leaves = parent.Children;
            dist = zeros(size(leaves,1),1);
            for h = 1 : size(leaves,1)
                dist(h) = leaves(h).RefVector.distance(point);
            end
            [~,idx] = sort(dist,'ascend');
            winner1 = leaves(idx(1));
            winner2 = leaves(idx(2));

            if epoch == 0 && flag == 0
                if parent.VoronoiCard > 100
                    idx = randperm(parent.VoronoiCard);
                    idx = idx(1:int64(parent.VoronoiCard/10));
                    matrix = parent.voronoiMatrix();
                    d = pdist(matrix(idx,:));
                else
                    matrix = parent.voronoiMatrix();
                    d = pdist(matrix);
                end
                avg = mean(d);
                winner1.T = avg;
                winner2.T = winner1.T;
                parent.T = winner1.T;
                winner1.add2Voronoi(point);
                winner1.Time = winner1.Time + 1;
                winner1.RefVector.updateCoordW(point,EpsilonW);
                flag = 1;

            else
                if winner1.NeighboursSize > 1
                    % Use convex hull
                    xi = point.coordinates;
                    w1 = winner1.RefVector.coordinates;
                    wN = zeros(winner1.NeighboursSize, size(winner1.RefVector.coordinates,2));
                    for j = 1 : winner1.NeighboursSize
                        wN(j,:) = winner1.Neighbours(j).RefVector.coordinates;
                    end
                    spiegabile = checkDataIntoSimplex(xi, w1, wN) || winner1.RefVector.distance(point) < winner1.T;
                else
                    % Use threshold
                    spiegabile = winner1.RefVector.distance(point) < winner1.T;
                end
                if spiegabile
                % Dato spiegabile
                    winner1.add2Voronoi(point);
                    winner1.Time = winner1.Time + 1;
                    winner1.RefVector.updateCoordW(point,EpsilonW);
                    for j = 1 : winner1.NeighboursSize
                        neighbour = winner1.Neighbours(j);
                        neighbour.RefVector.updateCoordN(EpsilonN,point, EpsilonW, winner1);
                    end
                    winner1.addEdge(winner2);
                    %update of the threshold 
                    winner1.updateThreshold();
                    winner2.updateThreshold();
%                     if winner1.RefVector.distance(winner2.RefVector) > winner1.T
%                         winner1.T = winner1.RefVector.distance(winner2.RefVector);
%                     end
                else
                % Dato non spiegabile
                    node = Node(parent, point.coordinates, parent.Root);
                    node.add2Voronoi(point);
                    node.T = parent.T;
                    counter = 0;
                end
            end
        end 
    end
    
end