function [ f ] = plotCluster( M, leaves, pcaFlag, coordinates )
    f = figure('units','normalized');%,'outerposition',[0 0 1 1]);
    ax = axes();
    ax.XLim = [-5,5];
    ax.YLim = [-5,5];
    if pcaFlag == 1
        [coeff, score, ~] = pca(M);
        scatter(score(:,1), score(:,2), '.', 'r');
    else
        if size(M, 2) > 2
            plot3(M(:,1),M(:,2),M(:,3),'*','Color',[0.9290, 0.6940, 0.1250])
        else
            plot(M(:,1),M(:,2),'*','Color',[0.9290, 0.6940, 0.1250])
        end
    end
    hold on
    for h = 1 : size(leaves,2)
        if pcaFlag == 1
            coord = leaves(h).RefVector.coordinates * coeff;
        else
            coord = leaves(h).RefVector.coordinates;
        end
        if size(M, 2) > 2
            plot3(coord(1),coord(2),coord(3), 'o','LineWidth',1,...
                'MarkerFaceColor',[0, 0.4470, 0.7410],...
                'MarkerEdgeColor','none','MarkerSize',10)
        else
            plot(coord(1),coord(2), 'o','LineWidth',1,...
                'MarkerFaceColor',[0, 0.4470, 0.7410],...
                'MarkerEdgeColor','none','MarkerSize',10)
        end
%         text(coord(1), coord(2), num2str(leaves(h).VoronoiCard));
        if leaves(h).Parent ~= false
            for j = 1 : leaves(h).EdgesSize
                edge = leaves(h).Edges(j);
                if edge.node1.IsLeaf && edge.node2.IsLeaf
                    if pcaFlag == 1
                        coord1 = edge.node1.RefVector.coordinates*coeff;
                        coord2 = edge.node2.RefVector.coordinates*coeff;
                    else
                        coord1 = edge.node1.RefVector.coordinates;
                        coord2 = edge.node2.RefVector.coordinates;
                    end
                    if size(M, 2) > 2
                        plot3([coord1(1) coord2(1)], [coord1(2) coord2(2)], [coord1(3) coord2(3)],...
                            'Color', [0, 0.4470, 0.7410]	, 'LineWidth', 2 );
                    else
                        plot([coord1(1) coord2(1)], [coord1(2) coord2(2)],...
                            'Color', [0, 0.4470, 0.7410]	, 'LineWidth', 2 );
                    end
%                     text((coord1(1) + coord2(1))/2, (coord1(2) + coord2(2))/2, num2str(leaves(h).Edges(j).age));
                end
            end
        end
    end
    if nargin == 4
        pointcoord = coordinates * coeff;
        scatter(pointcoord(1), pointcoord(2), 50,'k','filled')
    end
    hold off

end

