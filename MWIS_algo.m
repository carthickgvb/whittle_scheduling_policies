%Implementation of Max Weight independent set 
function out = MWIS_algo(n,G,M,psi,l,C,sim_time,inputs)

X = zeros(n,1);
A = adjacency(G);
idx = zeros(1,10);

Cost = zeros(sim_time,1);
Transmissions = zeros(sim_time,1);
Packets_dropped = zeros(sim_time,1);

for slot = 1:sim_time
    idx = MWIS(n,A,X);
     
    Cost(slot) = 0;
    for m=1:n
        Cost(slot) = Cost(slot) + C*X(m) + idx(m)*f(psi(m)); 
    end
    
    Transmissions(slot) = 0;
    for m=1:n
        k = inputs(slot,m);
        Packets_dropped(slot) = Packets_dropped(slot) + max(0,X(m) - idx(m)*psi(m) + k - M(m)); 
        X(m) = min(X(m) - idx(m)*psi(m) + k,M(m));
        Transmissions(slot) = Transmissions(slot) + idx(m)*psi(m);
    end
    
end

out = struct('cost',Cost,'transmissions',Transmissions,'packets_dropped',Packets_dropped);


end