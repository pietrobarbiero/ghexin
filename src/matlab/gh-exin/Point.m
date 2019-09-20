classdef Point < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        coordinates
        node
        position
        number
        outlier
        %pdist
    end
    
    methods
        function point = Point(coord, number)
            if nargin >= 1 
                %point.pdist = 1;
                point.coordinates = coord;
                point.outlier = 0;
            end
            if nargin == 2
                point.number = number;
            end
        end
        
        function d = distance(x, y)
            % distance calculated according to the norm "p" (=pdist)
            %d = (sum(abs(x.coordinates - y.coordinates).^x.pdist)).^(1/x.pdist);
            % Euclidean distance
            d = sqrt(sum((x.coordinates - y.coordinates).^2));
            % CityBlock distance
            % d = sum(abs(x.coordinates - y.coordinates));
            % distance calculated according to the mahalanobis distance
            %d = mahal(x.coordinates,y.coordinates);
        end
        function updateCoordW(x, point, Epsilon)
            Deltacoord = Epsilon * (point.coordinates - x.coordinates)/x.node.Time; 
            x.coordinates = x.coordinates + Deltacoord;
        end
        
        function updateCoordN(x, EpsilonN, point, EpsilonW, winner1)
            if EpsilonN < 0
               % Epsilon N auto 
               %Perche voronoi card e non time???
               sigma = - EpsilonN/x.node.VoronoiCard;
               alfaN = EpsilonW * exp(-(winner1.RefVector.distance(x)^2)/(2*(sigma)));         
               Deltacoord = alfaN * (point.coordinates - x.coordinates)/winner1.Time; 
               x.coordinates = x.coordinates + Deltacoord;
            else 
               Deltacoord = EpsilonN * (point.coordinates - x.coordinates)/x.node.Time; 
               x.coordinates = x.coordinates + Deltacoord; 
            end
        end
        
        function added2Node(point, node, pos)
            point.node = node;
            point.position = pos;
        end
        
        function deleteNode(point)
            if point.node ~= 0
                point.node.deleteFromVoronoi(point, point.position);
                point.node = 0;
                point.position = 0;
            end
        end
    end
    
end

