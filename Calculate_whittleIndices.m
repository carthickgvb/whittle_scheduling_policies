function out = Calculate_whittleIndices(M,PSI,C,l,mis,GoC)
n = length(M);
w_idx = zeros(n,max(M)+1);
for m=1:n
    if GoC == "graph"
        output = EquationSolve_graph(M(m),PSI(m),C,l(m),mis(m));
    elseif GoC == "clique"
        output = EquationSolve_clique(M(m),PSI(m),C,l(m));
    else
        error("wrong argument");
    end
    for p=1:(M(m)+1)
        w_idx(m,p) = output(p);
    end
end

%Interpolate the calculated whittle index from the point it starts to
%decrease
increase_idx = zeros(1,n);
for i=1:n
    for j=2:M(i)+1
        if (w_idx(i,j) - w_idx(i,j-1) > 0)
            increase_idx(i) = j-1;
            w_idx(i,1:M(i)+1) = interp1(1:increase_idx(i),w_idx(i,1:increase_idx(i)),1:M(i)+1,'linear','extrap');
            break;
        end
    end
end
out = w_idx;
end