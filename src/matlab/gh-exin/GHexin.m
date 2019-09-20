function [ nodes, leaves, outliers] = GHexin(M, HMax, Hpercentage, EpsilonW, EpsilonN, AgeMax, profile)
    % DGSOT Summary
    % Detailed explanation goes here
    
    fprintf('Starting GHexin\n');  
    
    rowSize = size(M,1);
    
    % Initializing the points of the matrix
    points(rowSize) = Point();
    for i = 1 : rowSize
        points(i) = Point(M(i,:), i);
    end
    
    fprintf('Initializing the root of the tree\n');
    % - finding the centroid of the dataset
    avg = mean(M,1);
    root = Node(false, avg, 0);
    for i = 1 : rowSize
        root.add2Voronoi(points(i));
    end
    
    k = 0;    
    outliers = [];
    flag = 0;
    flag2 = 0;
    
    %%Entire cycle
    while(1)

%       Plot structure of the tree
%         if k > 0
%             plotGraph([root root.subTree()], root.SubTreeCard, 1);
%         end
        
        %%Leaves initialization
        k = k + 1;
        [Kleaves, newOutliers] = root.leavesKH(k, HMax, profile);
        outliers = [outliers; newOutliers];

        % Plot current clusters
%         plotCluster(M, Kleaves);
        
        %%Wait two times without children before exiting
        leavesSize = size(Kleaves,2);
        if leavesSize == 0
            if flag2 == 1
                break;
            else 
                flag2 = 1;
            end
        else 
            flag2 = 0;
        end 
        fprintf('Vertical Growth - height = %s\n', num2str(k));
        
        for i = 1 : leavesSize
            parent = Kleaves(i);
            child1 = Node(parent, parent.RefVector.coordinates, root);
            child2 = Node(parent, parent.RefVector.coordinates, root);
            
            counter = 0;
%             prevChildren = size(parent.Children,1);
            epoch = 0;
            avgHcc = realmax;
            while (flag == 0 || epoch < 2 || counter < 2 && avgHcc > Hpercentage * parent.Hcc())...
                     && size(parent.Children, 1) >= 2
%           while (counter == 0 && (epoch < 2 || avgHcc > Hpercentage * parent.Hcc()))...
%                    && size(parent.Children, 1) >= 2
                learningProcess(parent, EpsilonW, EpsilonN, epoch);
                deletedE = pruningEdge(parent, AgeMax);
                [pOutliers,deletedN] = pruningNode(parent, card);
                outliers = [outliers; pruningData(pOutliers, parent, root)];
                epoch = epoch + 1;
                flag = 1;
                avgHcc = parent.HChildren();
                fprintf('%d/%d leaves - epoch = %d - avgHcc = %d - children = %d - deleted node = %d - deleted edge = %d\n', i, leavesSize, epoch, avgHcc, size(parent.Children,1), deletedN, deletedE);
%               if size(parent.Children,1) > prevChildren
                if epoch == 1
                    prevAvg = parent.Hcc();
                end
                if avgHcc < prevAvg
%                     prevChildren = size(parent.Children,1);
                    prevAvg = avgHcc;
                else 
                    counter = counter + 1;
                end
            end
            
             % connected components
              connectedComponents(parent);
           if size(parent.Children) == 0                                                                                                                      
              parent.IsLeaf = true;
           end                                                                 
        end
    end
    nodes = [root root.SubTree];
    leaves = root.subtreeLeaves();
    
end

        