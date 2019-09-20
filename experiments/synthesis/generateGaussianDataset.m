function [ Y , compIdx] = generateGaussianDataset()
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here
    n = 500;
    
    c1 = 3;
    mu1 = [0.45, 0.45;...
            0.7, 0.6;...
            0.5, 0.7];
    sigma1 = cat(3, 0.005*eye(2,2),0.005*eye(2,2),0.005*eye(2,2));
    p1 = ones(1,c1)./c1;

    c2 = 4;
    mu2 = [1.3, 1.3;...
           1.1, 1.1;...
           1.3, 1.1;...
           1.1, 1.3];
    sigma2 = cat(3, 0.005*eye(2,2),0.005*eye(2,2),0.005*eye(2,2),0.005*eye(2,2));
    p2 = ones(1,c2)./c2;

    gm1 = gmdistribution(mu1,sigma1,p1);
    gm2 = gmdistribution(mu2,sigma2,p2);


    [Y1,compIdx1] = random(gm1,n*3);
    [Y2,compIdx2] = random(gm2,n*4);

    compIdx2 = compIdx2 + c1;

    Y = [Y1;Y2];
    compIdx = [compIdx1; compIdx2];
    
    figure
    scatter(Y(:,1),Y(:,2),10,'.') % Scatter plot with points of size 10
    hold on
    gmPDF1 = @(x,y)reshape(pdf(gm1,[x(:) y(:)]),size(x));
    fcontour(gmPDF1,'LineWidth',2)
    hold on
    gmPDF2 = @(x,y)reshape(pdf(gm2,[x(:) y(:)]),size(x));
    fcontour(gmPDF2,'LineWidth',2)
    hold off

end

