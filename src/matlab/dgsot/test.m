clc;
dbstop if error;
% 
% M = [1 0 0 2 2 1 13 12 22 7.6; 0 2 4 3 1 2 5 10 33 5.5; 2 3 4 3 1 22 8 5 5 2.5];
load('microarray.mat');
M = microarray(:, 1 : 203)';

load('crcResponsivenessTernary.mat');
classes = crcResponsivenessTernary(1:203) + 2;
alfa = 0.1;
sigma0 = 2;
hubHeterogenity = 100;
profileHeterogenity = 20;
epsilonAD = 0.05;
epsilonET = 0.01;
kappa = 1;
hetType = 1;
if hetType == 2
    purity = 0.6;
    efficiency = 0.01;
    hubHeterogenity = purity;
    profileHeterogenity = efficiency;
    M = M';
end


[nodes, leaves, outliers] = DGSOT(M, alfa, sigma0, hubHeterogenity, profileHeterogenity, epsilonAD, epsilonET, kappa, hetType, classes); 
nodesCard = size(nodes,2);

plotGraph(nodes, nodes(1).VoronoiCard);

% for i = 1 : size(leaves,2)
%     if leaves(i).VoronoiCard > 10
%         VoronoiMatrix = leaves(i).voronoiMatrix();
%         [tissuesNodes, tissuesLeaves] = DGSOT(VoronoiMatrix', 0.1, 2, 10000000, 30, 0.6, 0.5, 1);
%         figure = plotGraph(tissuesNodes, 203, crcResponsivenessTernary(1:203));
%     end
% end

save("nodes.mat", "nodes");
save("leaves.mat", "leaves");
save("outliers.mat", "outliers");


