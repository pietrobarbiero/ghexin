function [ nodes, leaves, outliers, points] = DGSOT(M1, alfa, sigma0, distortionHeterogenity, profileHeterogenity, epsilonAD, epsilonET, kappa, hetType, classes, MaxLevels)
    % DGSOT Summary of this function goes here
    % Detailed explanation goes here
    
    outliers = [];
    columnsSize = size(M1,2);
    
    % Initializing the points of the matrix
    points(columnsSize) = Point_dgsot();
    for i = 1 : columnsSize
        points(i) = Point_dgsot(M1(:,i), i);
    end
    
    M = M1; 
    columns = size(M, 2);
    pointsNotOutlier = points;
    
    % Initializing the root of the tree -
    % - finding the centroid of the dataset
    avg = mean(M,2);
    root = Node_dgsot(false, avg,0);
    
    %Calculating the mean distance
    D = fractional_norm(M,1);
    meanDist1 = mean(D);
    meanDist = mean(meanDist1);
    
    % Assigning all points to the root
    for i=1:columns
        root.add2Voronoi(pointsNotOutlier(i));
    end
    
    % Initialization of variables
    vertFlag = true;
    height = 0;
    while vertFlag
        
        % Plot current clusters
        % There may be some nodes with no edges due to connected components
%         plotCluster(M, root.subtreeLeaves(), 0);

        MaxLevels = MaxLevels - 1;
        if MaxLevels < 0
            break;
        end
        
        % Vertical Growth
        height=height+1;
        fprintf('Vertical Growth - height = %s\n', num2str(height));
        Kleaves = root.leavesK(height);
        leavesSize = size(Kleaves,2);
        if leavesSize==0
            break;
        end 
        
        vertFlag = false;
        for i=1:leavesSize
            % Checking if the leaf has no Point_dgsot assinged
            if Kleaves(i).VoronoiCard == 0
                continue;
            % Delete child
            % Kleaves(i).Parent.deleteChild(Kleaves(i));
            % Verifying heterogenity of the leaf i
            elseif Kleaves(i).distortion() > distortionHeterogenity && Kleaves(i).VoronoiCard > profileHeterogenity
                vertFlag = true;
                parent = Kleaves(i);
                
                % Creation of 2 children
                Node_dgsot(parent, parent.RefVector.coordinates, root);
                Node_dgsot(parent, parent.RefVector.coordinates, root);
               
                % Learning Process associated to the vertical growth 
                learningProcess(parent, alfa, sigma0, epsilonET, kappa, meanDist);
                
                % Horizontal Growth
%                 fprintf('Horizontal Growth\n');
                horizontalGrowth(parent, epsilonAD, epsilonET, sigma0, alfa, kappa, meanDist);
            end     
        end
    end
    
%   Plot final clusters    
%     figure
%     scatter3(M(1,:),M(2,:), M(3,:), '.','r')
%     hold on
%     leaves = root.subtreeLeaves();
%     for i=1:size(leaves,2)
%         coord=leaves(i).RefVector.coordinates;
%         scatter3(coord(1),coord(2),coord(3), '*','y')
%         text(coord(1), coord(2), coord(3), num2str(leaves(i).VoronoiCard));
%     end
%     hold off

    nodes = [root, root.subTree()];
    leaves = root.subtreeLeaves();
end

        