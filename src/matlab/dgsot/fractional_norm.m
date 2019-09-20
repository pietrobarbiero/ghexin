function D = fractional_norm(x1,x2,p)

    if nargin == 3
        D = sum(abs(x1-x2).^(p))^(1/p);
    elseif nargin == 2
        M = x1;
        p = x2;
        D = zeros(size(M,2));
        parfor i = 1 : size(M,2)
            v = zeros(size(M,2),1);
            for j = 1 : size(M,2)
                v(j) = sum(abs(M(:,i)-M(:,j)).^(p))^(1/p);
            end
            D(:,i) = v';
        end
    end