function [ nodes, leaves, outliers, points, stats] = GHexinOldHcc(M, HMax, Hpercentage, EpsilonW, EpsilonN, AgeMax, card, maxHeight, plot, min_epochs, avgT)
    % Spiegazione dei parametri
    % M - è la matrice su cui si effettua il clustering gerarchico
    % HMAX - è il parametro di stop. E' il valore massimo  di HCC che i cluster
    % possono avere. Provato: 0.001 - 0.01 - 0.001: molto restrittivo;uso 0.4
    % Hpercentage - è il parametro che regola la percentuale di miglioramento
    % di HCC per livello. Più è basso più epoche sono effettuate, più figli
    % generati ad ogni livellol padre. Provato 0.9: fa
    % fare pochi figli; 0.5 - 0.7: buono ma non esce mai per questo; 0.8 meglio?
    % EpsilonW -  gestisce la motilità e la velocità di apprendimento dei nodi.
    % Inoltre controlla anche la quantità di nodi che vengono creati ad ogni
    % iterazione
    % Valori provati: 0.1 motilità insufficiente; 0.8 troppo forse; 0.5 buona
    % EpsilonN - come sopra. Valori provati 0.01: insufficiente; 0.1 forse troppo
    % 0.05: buono. Valore negativo -> valore data driven NON FUNZIONA
    % AgeMax - età massima degli edge. All'aumentare dell'età si creano più edge e meno nodi
    % valore che regola in parte il numero di outliers. Valore provato:
    % 10: restrittivo. Si potrebbe condizionare al numero di punti.
    % Meglio: 5 (M3), 3 (M1 - M2). Agemax = -1 -> valore data driven - NON FUNZIONA 
    % card - indica il numero di punti di un cluester sotto cui si smette di
    % dividere. Non da garanzie sulla numerosità ma valori più alti
    % garantiscono comunque una cardinalità maggiore per i cluster. A discapito
    % di HCC. Valori provati: 3 - si potrebbe alzare. 5 usato
    % MaxHeigth: Numero massimo di livelli di discesa
    
    stats = [];
    
    fprintf('Starting GHexin\n');  
    k = 1;    
    outliers = [];
    flag2 = 0; 
    
    % Initializing the points of the matrix
    rowSize = size(M,1);
    points(rowSize) = Point();
    for i = 1 : rowSize
        points(i) = Point(M(i,:), i);
    end  
    
    % Finding the centroid of the dataset
    avg = mean(M,1);
    root = NodeOldHcc(false, avg, 0);       
    for i = 1 : rowSize
        root.add2Voronoi(points(i));
    end 
    
    % Entire cycle
    while(k < maxHeight)
%       % Plot structure of the tree - DEBUGGING
%         if k > 0
%             plotGraph([root root.subTree()], root.SubTreeCard, 1);
%         end

        % Leaves initialization
        [Kleaves, newOutliers] = root.leavesKH(k, HMax, card);
        outliers = [outliers; newOutliers];
        k = k + 1;
       
        % Plot current clusters
        % There may be some nodes with no edges due to connected components
        if plot == 1
            plotCluster(M, root.subtreeLeaves(), 0); 
        end
        
        % Wait two times without children before exiting
        leavesSize = size(Kleaves,2);
        if leavesSize == 0
            if flag2 == 1
                break;
            else 
                flag2 = 1;
            end
        else 
            flag2 = 0;
        end 
        
        fprintf('Vertical Growth - height = %s\n', num2str(k-1));
        for i = 1 : leavesSize
            parent = Kleaves(i);
            child1 = NodeOldHcc(parent, parent.RefVector.coordinates, root);
            child2 = NodeOldHcc(parent, parent.RefVector.coordinates, root);
            counter = 0;
            epoch = 0;
            avgHcc = realmax;
            while ((epoch < min_epochs || avgHcc > Hpercentage * parent.Hcc()))...
                   && size(parent.Children, 1 ) >= 2 && counter < 2
                [~, stats_t] = learningProcessOldHcc(parent, EpsilonW, EpsilonN, epoch, AgeMax, avgT);
                
                stats = [stats stats_t];
                
%                 deletedE = pruningEdge(parent, AgeMax);
%                 [pOutliers,deletedN] = pruningNode(parent, card);
                [newOutliers, deletedN] = pruningNode(parent, card);
%                 newOutliers = pruningData(pOutliers, parent, root);
                outliers = [outliers; newOutliers];
                epoch = epoch + 1;
%                 if epoch > 1
                avgHcc = parent.HChildren();
%                 end
                fprintf('%d/%d leaves - epoch = %d - avgHcc = %d - children = %d - deleted node = %d - deleted points = %d\n', i, leavesSize, epoch, avgHcc, size(parent.Children,1), deletedN, size(newOutliers, 1));
                if epoch == 1
                    prevAvg = parent.Hcc();
                end
                if avgHcc ~= prevAvg
                    prevAvg = avgHcc;
                else 
                    counter = counter + 1;
                end
            end
            % Connected components
            connectedComponentsOldHcc(parent);
            if size(parent.Children) == 0
                parent.IsLeaf = true;
            end
        end
        %Show for each expanded leaf children generation - DEBUGGING
%         if plot == 1
%             plotCluster(parent.voronoiMatrix, parent.Children', 0);
%             for l = 1 : size(outliers,1)
%                 hold on
%                 coord=outliers(l).coordinates;
%                 scatter(coord(1),coord(2),'.','k');
%             end
%             hold off
%         end

    end
    nodes = [root root.SubTree];
    leaves = root.subtreeLeaves();
    
    if plot == 1
        plotCluster(M, root.subtreeLeaves(), 0); 
    end
end
