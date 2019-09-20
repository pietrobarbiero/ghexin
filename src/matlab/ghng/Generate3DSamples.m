function [Samples]=Generate3DSamples(DatasetNdx,SamplesNumber)
% Generate samples from some 3D datasets
% Inputs:
%   DatasetNdx=Index of the dataset (1 to 8)
%   SamplesNumber=Number of samples to generate
% Outpu:
%   Samples=Generated samples (one sample per column)
% Taken from Todd Wittman, MANIfold Learning Matlab Demo,
% http://www.math.ucla.edu/~wittman/mani/index.html

N=SamplesNumber;
ExParam=1;
switch DatasetNdx
    case 1  % Swiss Roll
        tt = (3*pi/2)*(1+2*rand(1,N));  
        height = 21*rand(1,N);
        Samples = [tt.*cos(tt); height; ExParam*tt.*sin(tt)]';
%         handles.ColorVector = tt';
%         updateString{1} = 'Swiss Roll example loaded.';
    case 2  % Swiss Hole
        % Swiss Roll with hole example taken from Donoho & Grimes
        tt = (3*pi/2)*(1+2*rand(1,2*N));  
        height = 21*rand(1,2*N);
        kl = repmat(0,1,2*N);
        for ii = 1:2*N
            if ( (tt(ii) > 9)&(tt(ii) < 12))
                if ((height(ii) > 9) & (height(ii) <14))
                    kl(ii) = 1;
                end;
            end;
        end;
        kkz = find(kl==0);
        tt = tt(kkz(1:N));
        height = height(kkz(1:N));
        Samples = [tt.*cos(tt); height; ExParam*tt.*sin(tt)]';     
%         handles.ColorVector = tt';
%         updateString{1} = 'Swiss Hole example loaded.';
    case 3  % Corner Planes
        k = 1;
        xMax = floor(sqrt(N));
        yMax = ceil(N/xMax);
        cornerPoint = floor(yMax/2);
        for x = 0:xMax
            for y = 0:yMax
                if y <= cornerPoint
                    X(k,:) = [x,y,0];
                    ColorVector(k) = y;
                else
                    X(k,:) = [x,cornerPoint+(y-cornerPoint)*cos(pi*ExParam/180),(y-cornerPoint)*sin(pi*ExParam/180)];
                    ColorVector(k) = y;
                end;
                k = k+1;
            end;
        end;
        Samples = X;
%         handles.ColorVector = ColorVector';
%         updateString{1} = 'Corner Planes example loaded.';
    case 4  % Punctured Sphere by Saul & Roweis
        inc = 9/sqrt(N);   %inc = 1/4;
        [xx,yy] = meshgrid(-5:inc:5);
        rr2 = xx(:).^2 + yy(:).^2;
        [tmp ii] = sort(rr2);
        Y = [xx(ii(1:N))'; yy(ii(1:N))'];
        a = 4./(4+sum(Y.^2));
        Samples = [a.*Y(1,:); a.*Y(2,:); ExParam*2*(1-a)]';
%         handles.ColorVector = handles.X(:,3);
%         updateString{1} = 'Punctured Sphere example loaded.';
    case 5  % Twin Peaks by Saul & Roweis
        inc = 1.5 / sqrt(N);  % inc = 0.1;
        [xx2,yy2] = meshgrid(-1:inc:1);
        zz2 = sin(pi*xx2).*tanh(3*yy2);
        xy = 1-2*rand(2,N);
        Samples = [xy; sin(pi*xy(1,:)).*tanh(3*xy(2,:))]';
        Samples(:,3) = ExParam * Samples(:,3);
%         handles.ColorVector = handles.X(:,3);
%         updateString{1} = 'Twin Peaks example loaded.';
    case 6  % 3D Clusters
%         numClusters = ExParam;
        numClusters = 3;
        numClusters = max(1,numClusters);
        Centers = 10*rand(numClusters,3);
        D = L2_distance(Centers',Centers',1);
        minDistance = min(D(find(D>0)));
        k = 1;
        N2 = N - (numClusters-1)*9;
        for i = 1:numClusters
            for j = 1:ceil(N2/numClusters)
               X(k,1:3) = Centers(i,1:3)+(rand(1,3)-0.5)*minDistance/sqrt(12);
               ColorVector(k) = i;
               k = k + 1;
           end;
           % Connect clusters with straight line.
           if i < numClusters
               for t = 0.1:0.1:0.9
                    X(k,1:3) = Centers(i,1:3) + (Centers(i+1,1:3)-Centers(i,1:3))*t;
                    ColorVector(k) = 0;
                    k = k+1;
                end;
           end;
        end;
        Samples = X;
%         handles.ColorVector = ColorVector;
%         updateString{1} = '3D Clusters example loaded.';
    case 7  % Toroidal Helix by Coifman & Lafon
        noiseSigma=0.05;   %noise parameter
        t = (1:N)'/N;
        t = t.^(ExParam)*2*pi;
        Samples = [(2+cos(8*t)).*cos(t) (2+cos(8*t)).*sin(t) sin(8*t)]+noiseSigma*randn(N,3);
%         handles.ColorVector = t;
%         updateString{1} = 'Toroidal Helix example loaded.';
    case 8  % Gaussian randomly sampled
        X = ExParam * randn(N,3);
        X(:,3) = 1 / (ExParam^2 * 2 * pi) * exp ( (-X(:,1).^2 - X(:,2).^2) / (2*ExParam^2) );
        Samples = X;
%         handles.ColorVector = X(:,3);
%         updateString{1} = 'Gaussian example loaded.';
    case 9  % Ueda's spiral
        epsilon = 0.5;
        t = (([0:N-1]/(N-1))*4*pi)';
        Samples = [(13-(0.5*t)).*cos(t) -(13-(0.5*t)).*sin(t) t] + epsilon*randn(N,3);        
end

Samples=Samples';

% --- L2_distance function
% Written by Roland Bunschoten, University of Amsterdam, 1999
function d = L2_distance(a,b,df)
if (size(a,1) == 1)
  a = [a; zeros(1,size(a,2))]; 
  b = [b; zeros(1,size(b,2))]; 
end
aa=sum(a.*a); bb=sum(b.*b); ab=a'*b; 
d = sqrt(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
d = real(d); 
if (df==1)
  d = d.*(1-eye(size(d)));
end

