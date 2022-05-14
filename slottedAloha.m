% %take input graph G and M and psi  
% n = 10;                                                                     %number of nodes
% 
% r = rand(n,2);
% 
% dist = zeros(n);
% 
% for i=1:n
%     for j=1:n
%         dist(i,j) = sqrt((r(i,1) - r(j,1))^2 + (r(i,2) - r(j,2))^2);
%         if(dist(i,j) > 0.3)                                                 %to vary the connectivity based on distance
%             dist(i,j) = 0;
%         end
%     end
% end
% 
% G = graph(dist);
% 
% p = plot(G,'Layout','force','EdgeLabel',G.Edges.Weight);
% 
% l = zeros(n,1);
% M = randi([2,6],n,1);
% psi = zeros(n,1);
% for i=1:n
%     psi(i) = randi([1,M(i)]);
%     l(i) = randi([1,M(i)]);
% end
% 
% %initialize
% X = zeros(n,1); %jobs in queue
% for m=1:n
%     X(m) = randi([0,M(m)]);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      inputs
function out = slottedAloha(n,G,X,M,psi,l,C,sim_time)


p =0.37;

input = zeros(sim_time,n);

Cost = zeros(sim_time,1);
Transmissions = zeros(sim_time,1);
Packets_dropped = zeros(sim_time,1);

for slot = 1:sim_time
    prob_id = zeros(n,1);       %captures the probability of tranmition
    active_id = zeros(n,1);     %denotes if node is ready to transmit
    transmit_id = zeros(n,1);   %denotes whether transmition occured or not
    
    for m=1:n
        active_id(m) = (X(m) >= psi(m));    %assuming psi as frame size
        if active_id(m)==1
            prob_id(m) = binornd(1,p);
        end
    end
    
    for m=1:n
        if prob_id(m) == 1
            collision = 0;
            nbr = neighbors(G,m);
            
            for j=1:size(nbr,1)
                if transmit_id(nbr(j)) == 1
                    collision = 1;
                    break;
                end
            end
            
            if collision == 0
                transmit_id(m) = 1;
            end
            
        end
    end
    
    Cost(slot) = 0;
    for m=1:n
        Cost(slot) = Cost(slot) + C*X(m) + prob_id(m)*f(psi(m)); %including cost of transmission during collision
    end
    
    Transmissions(slot) = 0;
    for m=1:n
        k = poissrnd(l(m));
        input(slot,m) = k; %store the arrivals at each slot
        Packets_dropped(slot) = Packets_dropped(slot) + max(0,X(m) - transmit_id(m)*psi(m) + k - M(m)); 
        X(m) = min(X(m) - transmit_id(m)*psi(m) + k,M(m));
        Transmissions(slot) = Transmissions(slot) + transmit_id(m)*psi(m);
    end
    
end

out = struct('cost',Cost,'transmissions',Transmissions,'inputs',input,'packets_dropped',Packets_dropped);


end