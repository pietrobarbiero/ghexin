classdef Edge < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        age
        ageMax
        node1
        node2
        mark
    end
    
    methods
        function edge = Edge(n1, n2)
            if nargin == 2
                edge.node1 = n1;
                edge.node2 = n2;
                edge.age = 0;
            end
            edge.mark = 0;
        end
        
        function deleteNode(edge)
            edge.node1.deleteEdge(edge);
            edge.node2.deleteEdge(edge);
            edge.delete();
        end
    end
    
end

