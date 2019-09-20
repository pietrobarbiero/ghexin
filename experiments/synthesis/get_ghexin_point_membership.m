function [ point_membership, min_class ] = get_ghexin_point_membership( points)
%GETDGSOTMEMBERSHIP Summary of this function goes here
%   Detailed explanation goes here
    n_points = numel(points);
    point_membership = zeros(n_points,1);
    for i = 1: n_points
        if points(i).outlier == 1
            point_membership(i) = NaN;
        else
            point_membership(i) = points(i).node.Number;
        end
    end
    min_class = min(point_membership)
    point_membership(:) = point_membership(:) - min(point_membership) + 1;

end

