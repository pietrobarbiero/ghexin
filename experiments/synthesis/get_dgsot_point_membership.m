function [  point_membership, min_class ] = get_dgsot_point_membership( points)
%GETDGSOTMEMBERSHIP Summary of this function goes here
%   Detailed explanation goes here
    n_points = numel(points);
    point_membership = zeros(n_points,1);
    for i = 1: n_points
        point_membership(i) = points(i).node.Number;
    end
    min_class = - min(point_membership);
    point_membership(:) = point_membership(:) - min_class + 1;
    
end

