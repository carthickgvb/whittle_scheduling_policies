function out = Active_ID(G,whittle_idx)

n = numnodes(G);
active_id = ones(n,1)*-1;

i= 0;
while(any(active_id == -1))
    i = mod(i,n) + 1;
    if (active_id(i) == -1)
        smallest = 1;
        nbr = neighbors(G,i);

        for j = 1:size(nbr,1)
            if (active_id(nbr(j)) == -1)
                if (whittle_idx(i) > whittle_idx(nbr(j)))
                    smallest = 0;
                    break
                end
            end
        end

        if (smallest == 1)
            active_id(i) = 1;
            for j = 1:size(nbr,1)
                active_id(nbr(j)) = 0;
            end
        end
    end
end

out = active_id;
end