function [ ] = connectedComponentsOldHcc( parent )
    nodes = parent.Children;
    nodesCard = size(nodes,1);  
    s = NaN(nodesCard,1)';
    t = NaN(nodesCard,1)';
    m = 1;
    for l = 1 : nodesCard
        nodes(l).ConnComp = l;
    end
    for l = 1 : nodesCard
        for n = 1 : nodes(l).EdgesSize
            if nodes(l).Edges(n).mark == 0
                s(m) = nodes(l).Edges(n).node2.ConnComp;
                t(m) = nodes(l).Edges(n).node1.ConnComp;
                m = m + 1;
                nodes(l).Edges(n).mark = 1;
            end
        end
    end
    s = s(~isnan(s));    
    t = t(~isnan(t));
    G = graph(s, t);
    ccomp = conncomp(G);

%     figure
%     plot(G, 'Layout', 'layered');

    cc = unique(ccomp);
    if size(cc,2) > 1
%         plotGraph([parent.Root parent.Root.SubTree], parent.VoronoiCard, 1);
        for l = 1 : size(cc,2)
            Coordsum = 0;
            VoronoiSum = 0;
            idx = find(ccomp==l);
            for m = 1 : size(ccomp,2)
                if ismember(m,idx)
                    Coordsum = Coordsum + nodes(m).RefVector.coordinates;
                    VoronoiSum = VoronoiSum + nodes(m).VoronoiCard;
                end
            end
            avg = Coordsum/size(ccomp==l,2);
            node = NodeOldHcc(parent, avg, parent.Root);
            node.T = parent.T;
            node.VoronoiCard = VoronoiSum;
            node.Highlight = 1;
            parent.Highlight = 1;
            for m = 1 : size(ccomp,2)
                if ismember(m,idx)
                    nodes(m).Parent.deleteChild(nodes(m));
                    node.addChild(nodes(m));
                end
            end
        end
%          plotGraph([parent.Root parent.Root.SubTree], parent.Root.VoronoiCard, 1);
    end
%     figure
%     plot(G, 'Layout', 'layered');

end

