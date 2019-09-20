function [  ] = learningProcess( parent, alfa, sigma0, epsilon, kappa, meanDist)
%     fprintf('Learning Process\n');

    % KLD: get all leaves that takes part in the learning process 
    if kappa >= parent.Height
        kappa = parent.Height-1;
    end
    parentk = parent;
    for k=1:kappa
        parentk=parent.Parent;
    end
    leaves = parentk.subtreeLeaves();

    % Cycle of epochs
    relativeError = epsilon+1;
    epoch = 0;
    maxEpoch=10000;
    error = zeros(maxEpoch,1);
    
    %   Relative error of first epoch is the parent error
    while relativeError>epsilon && epoch<maxEpoch        
        epoch= epoch+1;
        error(epoch)=0;
        
%         fprintf("Epoch: %d\n", epoch);
        % Voronoi Matrix --> covariance Sigma --> Mahalanobis distance
        %VoronoiMatrix = parent.voronoiMatrix();
        %Sigma = cov(VoronoiMatrix);
        %InvSigma = inv(Sigma);
        
        % Epoch
        for j = 1 : parent.VoronoiCard
            
            point = parent.VoronoiSet(j);
            dist = realmax;
            
            % Deleting prevoius point assignation
            point.deleteNode();
            
            % Winner search for point j
            for h = 1 : size(leaves,2)
                dist2 = leaves(h).RefVector.distance(point);
                % Using Mahalonobis
                %dist2 = ((point.coordinates-leaves(h).RefVector.coordinates)/Sigma)*(leaves(h).RefVector.coordinates-point.coordinates)';
                if dist2 < dist
                    dist = dist2;
                    index =  h;
                end  
            end
            error(epoch) = error(epoch) + dist;
            winner = leaves(index);
            winner.add2Voronoi(point);
            winner.Time = winner.Time +1;
            
            % Update Reference Vectors
            siblings = winner.Parent.Children;
            for h = 1: size(siblings)
                %siblings(h).Time = siblings(h).Time + 1;
                siblings(h).RefVector.updateCoord(siblings(h).Time, alfa, sigma0, winner.RefVector, point, meanDist);
            end
        end
        if epoch ~= 1
            relativeError = ((error(epoch-1) - error(epoch))/error(epoch-1));
        end
    end
end