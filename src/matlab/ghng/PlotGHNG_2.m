function [Handle]=PlotGHNG_2(Model,MyAxis,Level)
% Plot a GHNG model in 2D
% E.J. Palomo
% Inputs:
%   Model = GHNG model
%   MyAxis = Vector of four components to scale the plot

Handle = [];
hold on
PlotGHNGRec(Model,MyAxis,Level);

function PlotGHNGRec(Model,MyAxis,Level)

if ~isempty(Model) && (Level <= 4),
    
    Colors = [0.4 0.2 0;1 0 0;1 0.65 0;0 1 1];
%     MyColor = Colors(Level,:);
    MyColor = [0, 0.4470, 0.7410];
        
    % Draw the child maps
    if Level > 1
        for NdxNeuro=1:numel(Model.Child), 
            PlotGHNGRec(Model.Child{NdxNeuro},MyAxis,Level-1);
        end
    else
        
        % Draw the neurons
    %     plot(Model.Means(1,:),Model.Means(2,:),'or','LineWidth',1,'MarkerFaceColor',MyColor,'MarkerSize',15*(1/Level));
    
        if size(Model.Means, 1) > 2
            plot3(Model.Means(1,:),Model.Means(2,:),Model.Means(3,:),...
                'or','LineWidth',1,'MarkerFaceColor',MyColor,'MarkerSize',10,...
                'MarkerEdgeColor','none');
        else
            plot(Model.Means(1,:),Model.Means(2,:),...
                'or','LineWidth',1,'MarkerFaceColor',MyColor,'MarkerSize',10,...
                'MarkerEdgeColor','none');
        end
    %     plot(Model.Means(1,:)*MyAxis(2),Model.Means(2,:)*MyAxis(4),'or','LineWidth',1,'MarkerFaceColor',MyColor,'MarkerSize',15*(1/Level));    
        % Draw the vertical connections
        NumNeurons = size(Model.Means,2);
        for NdxUnit=1:NumNeurons,
            if isfinite(Model.Means(1,NdxUnit))
                NdxNeighbors = find(Model.Connections(NdxUnit,:));
                for NdxMyNeigh=1:numel(NdxNeighbors)
    %                 line([Model.Means(1,NdxUnit) Model.Means(1,NdxNeighbors(NdxMyNeigh))],...
    %                     [Model.Means(2,NdxUnit) Model.Means(2,NdxNeighbors(NdxMyNeigh))],'Color',MyColor);
                    
                    if size(Model.Means, 1) > 2
                        line([Model.Means(1,NdxUnit) Model.Means(1,NdxNeighbors(NdxMyNeigh))],...
                                [Model.Means(2,NdxUnit) Model.Means(2,NdxNeighbors(NdxMyNeigh))],...
                                [Model.Means(3,NdxUnit) Model.Means(3,NdxNeighbors(NdxMyNeigh))],'Color',MyColor);                
                    else
                        line([Model.Means(1,NdxUnit) Model.Means(1,NdxNeighbors(NdxMyNeigh))],...
                                [Model.Means(2,NdxUnit) Model.Means(2,NdxNeighbors(NdxMyNeigh))],'Color',MyColor);                
                    end
                    
    %                 line([Model.Means(1,NdxUnit)*MyAxis(2) Model.Means(1,NdxNeighbors(NdxMyNeigh))*MyAxis(2)],...
    %                     [Model.Means(2,NdxUnit)*MyAxis(4) Model.Means(2,NdxNeighbors(NdxMyNeigh))*MyAxis(4)],'Color',MyColor);
                end
            end
        end
    end
    
    
end