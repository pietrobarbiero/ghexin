function [  ] = horizontalGrowth( parent, epsilonAD, epsilonET, sigma0, alfa, kappa, meanDist)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    

    horizFlag = true;
    averageDistortions = zeros(1000,1);
    j = 1;
    % Average distortion with two children
    averageDistortions(j)= parent.averageDistortion();

    while horizFlag 
        % Creation of a sibling
        j = j+1;
        
        nodej = Node_dgsot(parent, parent.RefVector.coordinates, parent.Root);

        % Learning Process associated to the horizontal growth
        learningProcess(parent, epsilonET, sigma0, alfa, kappa, meanDist);
        % Calculate the average distortion with the new child
        averageDistortions(j)= parent.averageDistortion();
        if ((averageDistortions(j-1)-averageDistortions(j))/averageDistortions(j-1)<epsilonAD) || size(nodej.VoronoiSet,1)==0
            parent.deleteChild(nodej);
            learningProcess(parent, epsilonET, sigma0, alfa, kappa, meanDist);
            horizFlag = false;
        end
        if averageDistortions(j)==0
            horizFlag = false;
        end
    end
   
end

