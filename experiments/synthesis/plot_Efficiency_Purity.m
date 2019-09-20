function [Efficiencies, Purities] = plot_Efficiency_Purity(pointClasses, pointMembership, alg, min_class)   
    %% Calculating Efficiency and purity 
    classNumbers = unique(pointClasses(~isnan(pointClasses)));
    clusterNumbers = unique(pointMembership(~isnan(pointClasses)));
    numClasses = size(classNumbers,1);
    numCluster = size(clusterNumbers,1);
    numPoints = size(pointMembership,1);
    Efficiencies = zeros(numCluster,1);
    Purities = zeros(numCluster,1);
    pClass = zeros(numCluster,1);
    eClass = zeros(numCluster,1);
    
    countPerClass = zeros(numCluster,numClasses);
    
    for i = 1 : numPoints
        class = find(pointClasses(i)==classNumbers);
        cluster = find(clusterNumbers==pointMembership(i));
        countPerClass(cluster,class) = countPerClass(cluster,class) + 1;
    end

    numPointsPerCluster = zeros(numCluster,1);
    for i = 1: numCluster
        numPointsPerCluster(i) = sum(pointMembership==clusterNumbers(i));
    end
    numPointsPerClass = zeros(numClasses,1);
    for i = 1: numClasses
        numPointsPerClass(i) = sum(pointClasses==classNumbers(i));
    end
    
    for i = 1:numCluster
        % Calculating purity and efficiencies as the maximum purities and
        % efficiencies between classes 
        purities = countPerClass(i,:)./numPointsPerCluster(i);
        efficiencies = countPerClass(i,:)./numPointsPerClass(:)';
        [Purities(i), pClass(i)] = max(purities);
        pClass(i) = classNumbers(pClass(i));
        [Efficiencies(i), eClass(i)] = max(efficiencies);
        eClass(i) = classNumbers(eClass(i));
    end
    
    clusterNumbers = clusterNumbers + min_class - 1;
    %Plot leaves Purity
    f1 = figure;
    bar(clusterNumbers, Purities, 0.6);
    xlabel('Row Clusters');
    ylabel('Purity');
    str = 'Purities of '+alg+' leaves';
    title(str);
    for i = 1 : numCluster
        text(clusterNumbers(i), Purities(i)+0.05, num2str(pClass(i)));
    end
%     saveas(f1,char(alg+"Purity.eps"), 'eps');
    box off
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    set(f1,'PaperSize',[5.5 5]); %set the paper size to what you want
    print(f1,char(alg+'_purity'),'-dpdf') % then print it
    print(f1,char(alg+'_purity'),'-dpng') % then print it


    %Plot leaves Efficiency
    f2 = figure;
    bar(clusterNumbers, Efficiencies, 0.6);
    xlabel('Row Clusters');
    ylabel('Efficiency');
    axis([-inf inf 0 1.1])
    str = 'Efficiency of '+alg+' leaves';
    title(str);
    for i = 1 : numCluster
        text(clusterNumbers(i), Efficiencies(i)+0.05, num2str(pClass(i)));
    end
%     saveas(f2,char(alg+'Efficiency.eps'), 'eps'); 
    box off
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    set(f2,'PaperSize',[5.5 4.5]); %set the paper size to what you want
    print(f2,char(alg+'_efficiency'),'-dpdf') % then print it
    print(f2,char(alg+'_efficiency'),'-dpng') % then print it
end