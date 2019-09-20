function [Handle]=PlotGHNG(Model,MyAxis)
% Plot a GHNG model in 2D
% E.J. Palomo
% Inputs:
%   Model = GHNG model
%   MyAxis = Vector of four components to scale the plot

Handle = [];
hold on
PlotGHNGRec(Model,MyAxis,1);

function PlotGHNGRec(Model,MyAxis,Level)

if ~isempty(Model) && (Level <= 4),
    
    Colors = [0.4 0.2 0;1 0 0;1 0.65 0;0 1 1];
    MyColor = Colors(Level,:);
        
    % Draw the child maps
%     for NdxNeuro=1:numel(Model.Child), 
%         PlotGHNGRec(Model.Child{NdxNeuro},MyAxis,Level+1);
%     end
        
    % Draw the neurons
%     plot(Model.Means(1,:),Model.Means(2,:),'or','LineWidth',1,'MarkerFaceColor',MyColor,'MarkerSize',15*(1/Level));
    plot(Model.Means(1,:)*MyAxis(2)+MyAxis(1),Model.Means(2,:)*MyAxis(4)+MyAxis(3),'or','LineWidth',1,'MarkerFaceColor',MyColor,'MarkerSize',15*(1/Level));
%     plot(Model.Means(1,:)*MyAxis(2),Model.Means(2,:)*MyAxis(4),'or','LineWidth',1,'MarkerFaceColor',MyColor,'MarkerSize',15*(1/Level));    
    % Draw the vertical connections
    NumNeurons = size(Model.Means,2);
    for NdxUnit=1:NumNeurons,
        if isfinite(Model.Means(1,NdxUnit))
            NdxNeighbors = find(Model.Connections(NdxUnit,:));
            for NdxMyNeigh=1:numel(NdxNeighbors)
%                 line([Model.Means(1,NdxUnit) Model.Means(1,NdxNeighbors(NdxMyNeigh))],...
%                     [Model.Means(2,NdxUnit) Model.Means(2,NdxNeighbors(NdxMyNeigh))],'Color',MyColor);
                line([Model.Means(1,NdxUnit)*MyAxis(2)+MyAxis(1) Model.Means(1,NdxNeighbors(NdxMyNeigh))*MyAxis(2)+MyAxis(1)],...
                    [Model.Means(2,NdxUnit)*MyAxis(4)+MyAxis(3) Model.Means(2,NdxNeighbors(NdxMyNeigh))*MyAxis(4)+MyAxis(3)],'Color',MyColor);                
%                 line([Model.Means(1,NdxUnit)*MyAxis(2) Model.Means(1,NdxNeighbors(NdxMyNeigh))*MyAxis(2)],...
%                     [Model.Means(2,NdxUnit)*MyAxis(4) Model.Means(2,NdxNeighbors(NdxMyNeigh))*MyAxis(4)],'Color',MyColor);
            end
        end
    end
    
    
end