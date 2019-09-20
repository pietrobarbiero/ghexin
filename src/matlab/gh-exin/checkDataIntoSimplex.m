function [ result ] = checkDataIntoSimplex( x_i, firstWinner, firstWinnerNeighbors)
%checkDataIntoSimplex check if x_i is into simplex made by first winnner and
%its neighbors
%%param x_i row vector of the data
%%param firstWinner row vector of the first winner weight with regard to x_i
%%param firstWinnerNeighbors matrix of the neighbors of firstWinner. each
%%row contains a neighbor weight

numberOfNeighbors = size(firstWinnerNeighbors, 1);
if (numberOfNeighbors > 1 ) %%simplex needs at least two points other than firstWinner
    v_firstWinner = firstWinner - x_i;
    v_firstWinnerNeighbors = firstWinnerNeighbors - repmat(x_i, numberOfNeighbors, 1);
    a = v_firstWinner + sum(v_firstWinnerNeighbors);

    p_firstWinner = dot(v_firstWinner, a);
    lastSign = sign(p_firstWinner);
    result = true;
    for i = 1:numberOfNeighbors
        p_i = dot(v_firstWinnerNeighbors(i, :), a);
        if( sign(p_i) ~= 0 )
            if( lastSign == 0 || sign(p_i) == lastSign )
                lastSign  = sign(p_i);
                result = false;
            else
                result = true;
                return;
            end
        end
    end % end for
else
    result = NaN;
    return;
end
end