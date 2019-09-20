function [ M, outliers] = antiHubFilter( M1, nNeigh )
%   AntiHub Filter in order to identify and remove possible outliers
%   from the original dataset

    columns = size(M1,2);
    neighbours = Inf(columns, nNeigh, 2);
    outliers = ones(columns, 1);
    for i = 1 : columns
        for j = 1 : columns
            if j ~= i
                dist = sum(abs(M1(:,i) - M1(:,j)));
                maxdist = [0,0];
                for k = 1 : nNeigh
                    if maxdist(1) <= neighbours(i,k,2)
                        maxdist(1) = neighbours(i,k,2);
                        maxdist(2) = k;
                    end
                end
                if dist <= maxdist(1)
                    neighbours(i,maxdist(2),2) = dist;
                    neighbours(i,maxdist(2),1) = j;
                end
            end
        end
        for k = 1 : nNeigh
            outliers(neighbours(i,k,1)) = 0;
        end
    end
    
    M = M1(:,outliers == 1); 
end

