classdef Node_dgsot < handle
    %Node_dgsot Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Children
        GeneNames
        Height
        IsLeaf
        Net
        Number
        Parent
        Performance
        Ph3Bicluster
        Ph3Hmin
        Ph3K
        Ph3Tissues
        RefVector
        Root
        SubTree
        SubTreeCard
        TableResults
        Time
        VoronoiSet
        VoronoiCard
    end
    
    methods
        function Node_dgsot = Node_dgsot(parent,coord, root)
            if nargin ~= 0
                if isequal(parent,false)
                    Node_dgsot.Parent = false;
                    Node_dgsot.Height = 1;
                    Node_dgsot.Root = Node_dgsot;
                    Node_dgsot.Number = 1;
                else
                    Node_dgsot.Parent = parent;
                    Node_dgsot.Height = parent.Height+1;
                    parent.addChild(Node_dgsot);
                    Node_dgsot.Root = root;
                    Node_dgsot.Number = Node_dgsot.Root.SubTreeCard + 1;
                end
                point = Point_dgsot(coord);
                Node_dgsot.RefVector = point;
                Node_dgsot.VoronoiSet = [];
                Node_dgsot.VoronoiCard = 0;
                Node_dgsot.Children = [];
                Node_dgsot.SubTree = [];
                Node_dgsot.SubTreeCard = 0;
                Node_dgsot.IsLeaf = true;
                Node_dgsot.Time = 1;
            end
        end
        
        function Node_dgsot = addChild(Node_dgsot,child)
            if Node_dgsot.Children == 0
                Node_dgsot.Children = child;
            else
                Node_dgsot.Children = [Node_dgsot.Children; child];
            end
            Node_dgsot.IsLeaf = false;
            Node_dgsot.add2Subtree(child);
            parent = Node_dgsot.Parent;
            while parent~=false
                parent.add2Subtree(child);
                parent = parent.Parent;
            end
        end

        function Node_dgsot = deleteChild(Node_dgsot,child)
            % Node_dgsot.Children = Node_dgsot.Children(~isequal(Node_dgsot.Children,child));
            newChildren = [];
            for i = 1 : size(Node_dgsot.Children)
                if isequal(Node_dgsot.Children(i), child)
                    newChildren = [Node_dgsot.Children(1:i-1); Node_dgsot.Children(i+1:size(Node_dgsot.Children,1))];
                end
            end
            Node_dgsot.Children = newChildren;
            Node_dgsot.deleteFromSubTree(child);
            parent = Node_dgsot.Parent;
            while parent~=false
                parent.deleteFromSubTree(child);
                parent = parent.Parent;
            end
        end
                
        function Node_dgsot = add2Subtree(Node_dgsot, node2)
            if Node_dgsot.SubTreeCard == 0
                Node_dgsot.SubTree = node2;
                Node_dgsot.SubTree(10) = Node_dgsot;
                Node_dgsot.SubTreeCard = 1;
            else
                % Re-allocazione della matrice (dimensione doppia)
                if Node_dgsot.SubTreeCard == size(Node_dgsot.SubTree,2)
                    Node_dgsot.SubTree(Node_dgsot.SubTreeCard*2) = Node_dgsot;
                end
                % Incremento numero di celle piene e aggiungo il nuovo punto
                Node_dgsot.SubTreeCard = Node_dgsot.SubTreeCard + 1;
                Node_dgsot.SubTree(Node_dgsot.SubTreeCard) = node2;
            end
        end        
        
        function Node_dgsot = add2Voronoi(Node_dgsot, point) 
            %if isempty(Node_dgsot.VoronoiSet)
            if Node_dgsot.VoronoiCard == 0
                Node_dgsot.VoronoiSet = point;
                Node_dgsot.VoronoiSet(10) = Point_dgsot;
                Node_dgsot.VoronoiCard = 1;
            else
                % Re-allocazione della matrice (dimensione doppia)
                if (Node_dgsot.VoronoiCard == size(Node_dgsot.VoronoiSet,2))
                    Node_dgsot.VoronoiSet(Node_dgsot.VoronoiCard*2) = Point_dgsot;
                end
                % Incremento numero di celle piene e aggiungo il nuovo punto
                Node_dgsot.VoronoiCard = Node_dgsot.VoronoiCard + 1;
                Node_dgsot.VoronoiSet(Node_dgsot.VoronoiCard) = point;
            end
            point.added2Node(Node_dgsot, Node_dgsot.VoronoiCard);
        end
        
        function Node_dgsot = deleteFromSubTree(Node_dgsot, node2)
            for i = Node_dgsot.SubTreeCard : -1 : 1
                if isequal(Node_dgsot.SubTree(i), node2)
                    Node_dgsot.SubTree(i) = Node_dgsot.SubTree(Node_dgsot.SubTreeCard);
                    Node_dgsot.SubTree(Node_dgsot.SubTreeCard) = Node_dgsot;
                    Node_dgsot.SubTreeCard = Node_dgsot.SubTreeCard -1; 
                    break
                end
            end
        end
        
        function Node_dgsot = deleteFromVoronoi(Node_dgsot, point, position)
            if Node_dgsot.IsLeaf
                % Controllo se un halving � necessario
                if (Node_dgsot.VoronoiCard == size(Node_dgsot.VoronoiSet,2)/4)
                    Node_dgsot.VoronoiSet = Node_dgsot.VoronoiSet(1:Node_dgsot.VoronoiCard*2);
                end
                % Elimino il punto dal set di Voronoi
                if Node_dgsot.VoronoiSet(position) == point                
                    Node_dgsot.VoronoiSet(position) = Node_dgsot.VoronoiSet(Node_dgsot.VoronoiCard);
                    Node_dgsot.VoronoiSet(position).position = position;
                    Node_dgsot.VoronoiSet(Node_dgsot.VoronoiCard) = Point_dgsot();
                    Node_dgsot.VoronoiCard = Node_dgsot.VoronoiCard - 1;
                end
            end
        end
        
        function subTree = subTree(Node_dgsot)
            subTree = Node_dgsot.SubTree(1:Node_dgsot.SubTreeCard); 
        end
        
        function LeavesAtK = leavesK(Node_dgsot, k)
            if k == Node_dgsot.Height
                LeavesAtK = Node_dgsot;
            else
                LeavesAtK(Node_dgsot.SubTreeCard) = Node_dgsot;
                j = 0;
                for i = 1 : Node_dgsot.SubTreeCard
                    if Node_dgsot.SubTree(i).IsLeaf && Node_dgsot.SubTree(i).Height == k
                        j = j+1;
                        LeavesAtK(j) = Node_dgsot.SubTree(i);
                    end
                end
                LeavesAtK = LeavesAtK(1:j);
            end
        end
        
        function leaves = subtreeLeaves(Node_dgsot)
            if Node_dgsot.IsLeaf
                leaves = Node_dgsot;
            else
                leaves(Node_dgsot.SubTreeCard) = Node_dgsot;
                k = 0;
                for i = 1 : Node_dgsot.SubTreeCard
                    if Node_dgsot.SubTree(i).IsLeaf
                       k = k + 1;
                       leaves(k) = Node_dgsot.SubTree(i);
                    end
                end
                leaves = leaves(1:k);
            end
        end
        
        function het = heterogenity(Node_dgsot, hub, profile, type, pointClasses)
            het = true;
            if type == 1
                if Node_dgsot.VoronoiCard < 2
                    het = false;
                else
                    %Calculating distances in the sub-matrix
                    D = fractional_norm(Node_dgsot.voronoiMatrix,1);

                    %Calculating hubness of the sub-matrix
%                     if size(D,1) > 5
%                         nk = 5;
%                     else 
%                         nk = size(D,1);
%                     end
%                     [~,~,Nk] = hubness(D,nk);
%                     index = sum(Nk(Nk>2*nk));
                    if Node_dgsot.VoronoiCard < profile % || index < hub
                        het = false;
                    end
                end
            elseif type == 2
                minPurity = hub;
                minEfficiency = profile;
                classes = zeros(1,3);
                for h = 1 : Node_dgsot.VoronoiCard
                   classes(pointClasses(Node_dgsot.VoronoiSet(h).number)) = classes(pointClasses(Node_dgsot.VoronoiSet(h).number)) + 1;
                end

                purities = classes./Node_dgsot.VoronoiCard; 
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
        
        function [distortion, num] = distortion(Node_dgsot)
            distortion = 0;
            num = Node_dgsot.VoronoiCard;
           if Node_dgsot.VoronoiSet ~= 0
               for i = 1: Node_dgsot.VoronoiCard
                    distortion = distortion + Node_dgsot.RefVector.distance(Node_dgsot.VoronoiSet(i));
               end
           end
        end
        
        function avgDistortion = averageDistortion(Node_dgsot)
            distortion = 0;
            num = 0;
            for i = 1 : size(Node_dgsot.Children)
                [iDistortion, iNum] = Node_dgsot.Children(i).distortion();
                distortion = distortion + iDistortion;
                num = num + iNum;
            end
            avgDistortion = distortion / num;
        end
        
        function VoronoiMatrix = voronoiMatrix(Node_dgsot)
             VoronoiMatrix = NaN(size(Node_dgsot.RefVector.coordinates,1),Node_dgsot.VoronoiCard);
             for j = 1: Node_dgsot.VoronoiCard
                VoronoiMatrix(:,j)= Node_dgsot.VoronoiSet(j).coordinates;
             end
        end
    end
end