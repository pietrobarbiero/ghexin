function [ counter, stats ] = learningProcessOldHcc(parent, EpsilonW, EpsilonN, epoch, ageMax, avgT )
    counter = 1;
    flag = 0;
%     ageMax = parent.VoronoiCard/ageMax;

    stats = zeros(3, 1);
    
    for i = 1 : parent.VoronoiCard
        assert(parent.VoronoiCard ~= 0)
        point = parent.VoronoiSet(i);
        
        % plotcluster - DEBUGGING
%         if rem (i,100) == 0
%             plotCluster(parent.voronoiMatrix, parent.Children', 1, point.coordinates);
%         end
        
        if point.outlier == 0
            % point.coordinates
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
            % First time subnet initialization
            if flag == 0
                flag = 1;
                matrix = parent.voronoiMatrix();
                d = pdist(matrix);
%                 if parent.Height == 1
%                     avg = log(mean(d));
%                 else
                    avg = mean(d);
%                 end
                if epoch == 0
                    winner1.T = avg;
                    winner2.T = winner1.T;
                    parent.T = winner1.T;
                end
                winner1.add2Voronoi(point);
                winner1.Time = winner1.Time + 1;
                winner1.RefVector.updateCoordW(point,EpsilonW);
                winner1.addEdge(winner2);
            else
                if winner1.NeighboursSize >= length(point.coordinates)
                    % Use convex hull
                    xi = point.coordinates;
                    w1 = winner1.RefVector.coordinates;
                    wN = zeros(winner1.NeighboursSize, size(winner1.RefVector.coordinates,2));
                    for j = 1 : winner1.NeighboursSize
                        wN(j,:) = winner1.Neighbours(j).RefVector.coordinates;
                    end
                    if avgT == true
                        spiegabile = checkDataIntoSimplex(xi, w1, wN) || winner1.RefVector.distance(point) < winner1.T;
                    else
                        spiegabile = checkDataIntoSimplex(xi, w1, wN);
                    end
                    
                    stats(1) = stats(1) + 1;
                    
                else
                    % Use threshold
                    spiegabile = winner1.RefVector.distance(point) < winner1.T;
                    
                    stats(2) = stats(2) + 1;
                    
                end
                if spiegabile
                % Dato spiegabile
                    winner1.add2Voronoi(point);
                    winner1.Time = winner1.Time + 1;
                    winner1.RefVector.updateCoordW(point,EpsilonW);
                    % Create or update edge between w1 and w2 and update all other winner edge ages
                    winner1.addEdge(winner2);
                    % Update of the threshold 
                    winner1.updateThreshold(avg);
%                     if winner1.RefVector.distance(winner2.RefVector) > winner1.T
%                         winner1.T = winner1.RefVector.distance(winner2.RefVector);
%                     end
                    winner1.pruningSingleNodeEdges(ageMax);
                    % Update neighbours coordinate and thresholds
                    for j = 1 : winner1.NeighboursSize
                        neighbour = winner1.Neighbours(j);
                        neighbour.RefVector.updateCoordN(EpsilonN,point, EpsilonW, winner1);
                        neighbour.updateThreshold(avg);
                    end
                else
                % Dato non spiegabile
                    node = NodeOldHcc(parent, point.coordinates, parent.Root);
                    node.add2Voronoi(point);
                    node.T = parent.T;
                    counter = 0;
                    
                    stats(3) = stats(3) + 1;
                end
            end
        end 
    end
    
end