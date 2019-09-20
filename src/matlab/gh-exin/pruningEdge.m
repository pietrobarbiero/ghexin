function [deleted] = pruningEdge( parent, AgeMax )
%PRUNINGEDGE: delete all the edges having an age above the maximum age
%   If the edge has been created too early and it has not been resetted
%   then it would probably means that once the 2 nodes were close to each
%   other but now they are no more --> they are not neighbours
    deleted = 0;
    flag = 0;
    if AgeMax == -1
        flag = 1;
        erroreTotale = parent.averageDistortion();
        meanAge = parent.meanAge();
    end
    for i = 1:size(parent.Children,1)
        child = parent.Children(i);
        j = 1;
        while j <= child.EdgesSize
            %%Automatically calculate agemax for each edge
            if flag == 1
                errore1 = child.Edges(j).node1.distortion();
                errore2 = child.Edges(j).node2.distortion();
                child.Edges(j).ageMax = 1/2*(errore1 + errore2)/erroreTotale*meanAge;
                AgeMax = child.Edges(j).ageMax;
            end
            if child.Edges(j).age > AgeMax
               child.Edges(j).deleteNode;
               deleted = deleted + 1;
            else
                j = j+1;
            end
        end
    end
end

