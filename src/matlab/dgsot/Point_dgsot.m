classdef Point_dgsot < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        coordinates
        node
        position
        number
        %pdist
    end
    
    methods
        function Point_dgsot = Point_dgsot(coord, number)
            if nargin >= 1 
                %Point_dgsot.pdist = 1;
                Point_dgsot.coordinates = coord;
            end
            if nargin == 2
                Point_dgsot.number = number;
            end
        end
        
        function d = distance(x, y)
            % distance calculated according to the norm "p" (=pdist)
            %d = (sum(abs(x.coordinates - y.coordinates).^x.pdist)).^(1/x.pdist);
            % Euclidean distance
            %d = sqrt(sum(x.coordinates - y.coordinates).^2));
            % CityBlock distance
            d = sum(abs(x.coordinates - y.coordinates));
            % distance calculated according to the mahalanobis distance
            %d = mahal(x.coordinates,y.coordinates);
        end
        
        function updateCoord(x, time, alfa, sigma0, y, Point_dgsot, meanDist)
            sigma = sigma0/time;
            dist = x.distance(y);
            eta = exp(-(dist/(2*sigma^2)));
            phi = alfa*eta;
            coorDist = Point_dgsot.coordinates-x.coordinates;
            Deltacoord = phi.*(coorDist);
            actualDist = fractional_norm(x.coordinates, x.coordinates + Deltacoord, 1);
            if meanDist/2 > actualDist
                x.coordinates = x.coordinates + Deltacoord;
            end
        end
        function added2Node(Point_dgsot, node, pos)
            Point_dgsot.node = node;
            Point_dgsot.position = pos;
        end
        
        function deleteNode(Point_dgsot)
            Point_dgsot.node.deleteFromVoronoi(Point_dgsot, Point_dgsot.position);
            Point_dgsot.node = 0;
            Point_dgsot.position = 0;
        end
    end
    
end

