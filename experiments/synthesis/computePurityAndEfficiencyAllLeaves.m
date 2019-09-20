function [rEfficiency, rPurity, eClass, pClass] =  computePurityAndEfficiencyAllLeaves(Model, labelVect, plotFlag)
%%Retrieve all the leaves of a tree and computes their purity and efficiency index w.r.t their labels
%param Model is the model yield by GHNG training
%param labelVect this vector contains the class label of each element of the training set
%param plotFlag is true plots the efficiency and purity

leaves = GetCentroidsGHNG(Model); %% find Tree leaves
winners = TestGHNG(leaves,Model.Samples); %% get Voronoi set for each leave as winner indices
 
indicesRange = 1:length(labelVect); %it used to convert logical values back into indices 
numberOfLeaves = size(leaves, 2);

rEfficiency = zeros(1, numberOfLeaves);
rPurity = zeros(1, numberOfLeaves);
eClass = zeros(1, numberOfLeaves);
pClass = zeros(1, numberOfLeaves);

for i = 1:numberOfLeaves
    voronoiSetIndices_i = indicesRange(winners == i);
%     rowClusteringVoronoiCard_i = length(rowClusteringVoronoiSetIndices_i);
        [rEfficiency(i), rPurity(i), eClass(i), pClass(i)] = computePurityAndEfficiencySingleNode(labelVect, voronoiSetIndices_i);
end


if (nargin>2 && plotFlag == true)
%% Plot leaves Purity and Efficiency
    alg = "GHNG";
% Plot leaves Purity
    f1 = figure;
    bar(rPurity, 0.6);
%     for i = 1 : size(firstLeaves)
%         text(i, rPurity(i), num2str(eClass(i)));
%     end
    xlabel('Row Clusters');
    ylabel('Purity');
    str = 'Purity of GHNG leaves';
    title(str);
%     saveas(f1,[experimentDir, fileMask, 'Purity'], 'fig'); 
%     saveas(f1,[experimentDir, fileMask, 'Purity'], 'png');
% 
    box off
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    set(f1,'PaperSize',[5.5 4]); %set the paper size to what you want
    print(f1,char(alg+'_purity'),'-dpdf') % then print it
    print(f1,char(alg+'_purity'),'-dpng') % then print it
    
%Plot leaves Efficiency
    f2 = figure;
    bar(rEfficiency, 0.6);
%     for i = 1 : size(firstLeaves)
%         text(i, rEfficiency(i), num2str(pClass(i)));
%     end
    xlabel('Row Clusters');
    ylabel('Efficiency');
    str = 'Efficiency of GHNG leaves';
    title(str);
    box off
    set(findall(gcf,'-property','FontSize'),'FontSize',14)
    set(f2,'PaperSize',[5.5 4]); %set the paper size to what you want
    print(f2,char(alg+'_efficiency'),'-dpdf') % then print it
    print(f2,char(alg+'_efficiency'),'-dpng') % then print it
end

end %end function