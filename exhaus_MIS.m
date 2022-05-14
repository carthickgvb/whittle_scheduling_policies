%Finds the maximum independent set in neighbours of node 'v' in graph 'G'
function out = exhaus_MIS(G,v)

nbr = neighbors(G,v);

G_temp = G;
idx = 1;
for i=1:numnodes(G)
    if (~ismember(i,nbr))
        ids(idx) = i;
        idx = idx + 1;
    end
end

%Create a subgraph of neighbors of node of interest
G_temp = rmnode(G_temp,ids);

%figure;
%plot(G_temp,'Layout','force','EdgeLabel',G_temp.Edges.Weight);

%out = findMIS(G_temp);
out = findMIS_new(G_temp);

end