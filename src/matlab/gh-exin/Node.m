classdef Node < handle
    %NODE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Children
        ConnComp
        Edges
        EdgesSize
        Neighbours
        NeighboursSize
        GeneNames
        Height
        Highlight
        HCC
        IsLeaf
        Number
        Parent
        Performance
        RefVector
        Root
        SubTree
        SubTreeCard
        TableResults
        T
        Time
        TissuesNames
        VoronoiSet
        VoronoiCard
    end
    
    methods
        function node = Node(parent, coord, root)
            if nargin ~= 0
                if isequal(parent,false)
                    node.Parent = false;
                    node.Height = 1;
                    node.Root = node;
                    node.Number = 1;
                    node.T = 0;
                else
                    node.Parent = parent;
                    node.Height = parent.Height+1;
                    parent.addChild(node);
                    node.Root = root;
                    node.Number = node.Root.SubTreeCard + 1;
                    node.T = parent.T;
                end
                point = Point(coord);
                node.RefVector = point;
                point.node = node;
                node.VoronoiSet = [];
                node.VoronoiCard = 0;
                node.Children = [];
                node.SubTree = [];
                node.SubTreeCard = 0;
                node.IsLeaf = true;
                node.Time = 1;
                node.NeighboursSize = 0;
                node.Neighbours = Node();
                node.EdgesSize = 0;
                node.Edges = Edge();
                node.HCC = 0;
                node.ConnComp = 0;
                node.Highlight = 0;
            end
        end
        
        function node = addChild(node,child)
            if node.Children == 0
                node.Children = child;
            else
                node.Children = [node.Children; child];
            end
            child.Parent = node;
            child.Height = node.Height + 1;
            node.IsLeaf = false;
            node.add2Subtree(child);
            parent = node.Parent;
            while parent ~= false
                parent.add2Subtree(child);
                parent = parent.Parent;
            end
            child.Height = node.Height + 1;
        end

        function node = deleteChild(node,child)
            newChildren = [];
            for i = 1 : size(node.Children)
                if isequal(node.Children(i), child)
                    if i == 1 && i == size(node.Children,1)
                        newChildren = [];
                    elseif i == 1
                        newChildren = node.Children(i+1:size(node.Children,1));
                    elseif i == size(node.Children,1)
                        newChildren = node.Children(1:i-1);
                    else
                        newChildren = [node.Children(1:i-1); node.Children(i+1:size(node.Children,1))];
                    end
                    flag = 1;
                    break;
                end
            end
            if flag == 1
                eSize = child.EdgesSize;
                for i = 1 : eSize
                    child.Edges(1).deleteNode();
                end
                node.Children = newChildren;
                node.deleteFromSubTree(child);
                parent = node.Parent;
                while parent ~= false
                    parent.deleteFromSubTree(child);
                    parent = parent.Parent;
                end
    %             It doesn't work with connected Components function
    %             child.delete();
            end
        end
                
        function node = add2Subtree(node, node2)
            if node.SubTreeCard == 0
                node.SubTree = node2;
                node.SubTree(10) = Node;
                node.SubTreeCard = 1;
            else
                % Re-allocazione della matrice (dimensione doppia)
                if node.SubTreeCard == size(node.SubTree,2)
                    node.SubTree(node.SubTreeCard*2) = Node;
                end
                % Incremento numero di celle piene e aggiungo il nuovo punto
                node.SubTreeCard = node.SubTreeCard + 1;
                node.SubTree(node.SubTreeCard) = node2;
            end
        end        
        
        function node = add2Voronoi(node, point) 
            %if isempty(node.VoronoiSet)
            if node.VoronoiCard == 0
                node.VoronoiSet = point;
                node.VoronoiSet(10) = Point;
                node.VoronoiCard = 1;
            else
                % Re-allocazione della matrice (dimensione doppia)
                if (node.VoronoiCard == size(node.VoronoiSet,2))
                    node.VoronoiSet(node.VoronoiCard*2) = Point;
                end
                % Incremento numero di celle piene e aggiungo il nuovo punto
                node.VoronoiCard = node.VoronoiCard + 1;
                node.VoronoiSet(node.VoronoiCard) = point;
            end
            point.added2Node(node, node.VoronoiCard);
        end
        
        function node = deleteFromSubTree(node, node2)
            for i = node.SubTreeCard : -1 : 1
                if isequal(node.SubTree(i), node2)
                    node.SubTree(i) = node.SubTree(node.SubTreeCard);
                    node.SubTree(node.SubTreeCard) = Node;
                    node.SubTreeCard = node.SubTreeCard -1; 
                    break
                end
            end
        end
        
        function node = deleteFromVoronoi(node, point, position)
            % Controllo se un halving � necessario
            if (node.VoronoiCard == size(node.VoronoiSet,2)/4)
                node.VoronoiSet = node.VoronoiSet(1:node.VoronoiCard*2);
            end
            % Elimino il punto dal set di Voronoi
            if node.VoronoiSet(position) == point
                if node.VoronoiCard > 1
                    node.VoronoiSet(position) = node.VoronoiSet(node.VoronoiCard);
                    if node.IsLeaf == 1
                        node.VoronoiSet(position).position = position;
                    end
                end
                node.VoronoiSet(node.VoronoiCard) = Point();
                node.VoronoiCard = node.VoronoiCard - 1;
            end
        end
        
        function subTree = subTree(node)
            subTree = node.SubTree(1:node.SubTreeCard); 
        end
        
        function LeavesAtK = leavesK(node, k)
            %Used to identify the root
            if k == node.Height
                LeavesAtK = node;
            else
                LeavesAtK(node.SubTreeCard) = Node;
                j = 0;
                for i = 1 : node.SubTreeCard
                    if node.SubTree(i).IsLeaf && node.SubTree(i).Height == k
                        j = j+1;
                        LeavesAtK(j) = node.SubTree(i);
                    end
                end
                LeavesAtK = LeavesAtK(1:j);
            end
        end
        
        function [LeavesAtKH, outliers]= leavesKH(root, k, H, profile)
            outliers = [];
            %Used to idetify the root
            if k == root.Height
                LeavesAtKH = root;
            else
                LeavesAtKH(root.SubTreeCard) = Node;
                j = 0;
                i = 0;
                while i < root.SubTreeCard
                    i = i + 1;
                    node = root.SubTree(i);
                    if node.IsLeaf && node.Height == k &&... 
                            node.Hcc() > H && node.VoronoiCard > profile
                        j = j+1;
                        LeavesAtKH(j) = node;
                    elseif node.VoronoiCard < 1
                        % It may be a parent which have not been able to
                        % create any child. Its point have been reassigned
                        % It is deleted
%                       fprintf("New outlier\n");
                        node.Parent.deleteChild(node);
                        outliers = [outliers; node.VoronoiSet(1:node.VoronoiCard)'];
                    end
                end
                LeavesAtKH = LeavesAtKH(1:j);
            end
        end
        
        function leaves = subtreeLeaves(node)
            if node.IsLeaf
                leaves = node;
            else
                leaves(node.SubTreeCard) = Node;
                k = 0;
                for i = 1 : node.SubTreeCard
                    if node.SubTree(i).IsLeaf
                       k = k + 1;
                       leaves(k) = node.SubTree(i);
                    end
                end
                leaves = leaves(1:k);
            end
        end
        
        function het = heterogenity(node, hub, profile, type, pointClasses)
            het = true;
            if type == 1
                if node.VoronoiCard < 2
                    het = false;
                else
                    %Calculating distances in the sub-matrix
                    D = fractional_norm(node.voronoiMatrix,1);

                    %Calculating hubness of the sub-matrix
                    if size(D,1) > 5
                        nk = 5;
                    else 
                        nk = size(D,1);
                    end
                    [~,~,Nk] = hubness(D,nk);
                    index = sum(Nk(Nk>2*nk));
                    if index < hub || node.VoronoiCard < profile
                        het = false;
                    end
                end
            elseif type == 2
                minPurity = hub;
                minEfficiency = profile;
                classes = zeros(1,3);
                for h = 1 : node.VoronoiCard
                   classes(pointClasses(node.VoronoiSet(h).number)) = classes(pointClasses(node.VoronoiSet(h).number)) + 1;
                end

                purities = classes./node.VoronoiCard; 
                % Calculating purity of the gene as the maximum purities
                % between classes
                purity = max(purities);
                class3 = count(num2str(pointClasses), '3');
                class2 = count(num2str(pointClasses), '2');
                class1 = count(num2str(pointClasses), '1');
                efficiencies(1) = classes(1)./class1;
                efficiencies(2) = classes(2)./class2;
                efficiencies(3) = classes(3)./class3;
                % Calculating purity of the gene as the maximum purities
                % between classes
                efficiency = max(efficiencies);
                if ((efficiency > minEfficiency)&& purity > minPurity)
                    het = false;
                end
            end
        end
        
        function [avg] = distortion(node)
            distortion = 0;
            if node.VoronoiCard ~= 0
               for i = 1: node.VoronoiCard
                    distortion = distortion + node.RefVector.distance(node.VoronoiSet(i));
               end
               avg = distortion / node.VoronoiCard; 
            else 
               avg = 0;
            end
        end
        
        function avgDistortion = averageDistortion(node)
            distortion = 0;
            num = 0;
            for i = 1 : size(node.Children)
%                 if node.Children(i).VoronoiCard > 1
                    iDistortion = node.Children(i).distortion();
                    distortion = distortion + iDistortion;
                    num = num + 1;
%                 end
            end
            avgDistortion = distortion / num;
        end
        
        function VoronoiMatrix = voronoiMatrix(node)
             VoronoiMatrix = NaN(node.VoronoiCard,size(node.RefVector.coordinates,2));
             for i = 1: node.VoronoiCard
                VoronoiMatrix(i,:)= node.VoronoiSet(i).coordinates;
             end
        end
        function hcc = Hcc(node)
            if node.HCC ~= 0 && node.IsLeaf == false
                hcc = node.HCC;
                return;
            end
            M = node.voronoiMatrix();
            nColumns = size(M,2);
            nRows = size(M,1);
            residue = zeros(nRows,nColumns);
            f = zeros(nRows,nColumns);
            avg = mean(mean(M));
            parfor i = 1 : nRows
                rowmean = mean(M(i,:));
                for j = 1 : nColumns
                    clmmean = mean(M(:,j));
                    residue(i,j) = M(i,j) - rowmean - clmmean + avg;
                    f(i,j) = M(i,j) - rowmean;
                end
            end
            hccVecchio = 1/(nRows*nColumns)*sum(sum(residue.^2));
            F = sum(sum(f.^2)/nColumns)/nRows;
            hcc = hccVecchio + 1/F;
            node.HCC = hcc;
            %fprintf('H2 index = %s\n', h2Cidx);
        end
        
        function avgH = HChildren(node)
            Hs = ones(size(node.Children,1),1);
            for i = 1 : size(node.Children,1) 
                Hs(i) = node.Children(i).Hcc();
            end
            avgH = mean(Hs);
        end
        
        function edge = addEdge(n1,n2)
            flag = 0;
            if n1.EdgesSize > size(n1.Parent.Children, 1)
                fprintf("Error\n");
            end
            for i = 1 : n1.EdgesSize
                if n1.Edges(i).node1 == n1 && n1.Edges(i).node2 == n2
                    edge = n1.Edges(i);
                    flag = 1;
                elseif n1.Edges(i).node2 == n1 && n1.Edges(i).node1 == n2
                    edge = n1.Edges(i);
                    flag = 1;
                else
                    n1.Edges(i).age = n1.Edges(i).age + 1;
                end
            end
            if flag
                edge.age = 0;
                return
            end
            %Create new edge
            n1.EdgesSize = n1.EdgesSize + 1;
            n2.EdgesSize = n2.EdgesSize + 1;
            edge = Edge(n1,n2);
            n1.Edges(n1.EdgesSize) = edge;
            n2.Edges(n2.EdgesSize) = edge;
            n1.NeighboursSize = n1.NeighboursSize + 1;
            n1.Neighbours(n1.NeighboursSize) = n2;
            n2.NeighboursSize = n2.NeighboursSize + 1;
            n2.Neighbours(n2.NeighboursSize) = n1;
        end
        
        function deleteEdge(node, edge)
            for i = 1 : node.EdgesSize
                if isequal(node.Edges(i), edge)
                    node.Neighbours = [node.Neighbours(1:i-1) node.Neighbours(i+1:node.NeighboursSize)];
                    node.NeighboursSize = node.NeighboursSize - 1;
                    node.Edges = [node.Edges(1:i-1) node.Edges(i+1:node.EdgesSize)];
                    node.EdgesSize = node.EdgesSize - 1;
                    break;
                end
            end
        end
        
        function avg = meanAge(node)
            ages = [];
            for i = 1: size(node.Children,1)
                newAges = [node.Children(i).Edges.age];
                ages = [ages; newAges'];
            end
            avg = mean(ages);
        end
        
        function updateThreshold(node, avgP)
            d = double(node.NeighboursSize);
            for i = 1 : node.NeighboursSize
                d(i) = node.RefVector.distance(node.Neighbours(i).RefVector);
            end
            avgN = mean(d);
            if avgN > avgP
                node.T = avgN;
            end
        end
        
        function pruningSingleNodeEdge( node, AgeMax )
        %delete all the edges of the winner 
            if AgeMax == -1
                flag = 1;
                erroreTotale = parent.averageDistortion();
                meanAge = parent.meanAge();
            end
            j = 1;
            while j <= node.EdgesSize
                %%Automatically calculate agemax for each edge
                if flag == 1
                    errore1 = node.Edges(j).node1.distortion();
                    errore2 = node.Edges(j).node2.distortion();
                    node.Edges(j).ageMax = 1/2*(errore1 + errore2)/erroreTotale*meanAge;
                    AgeMax = node.Edges(j).ageMax;
                end
                if node.Edges(j).age > AgeMax
                   node.Edges(j).deleteNode;
                else
                    j = j+1;
                end
            end
        end
    %%End of methods
    end
end