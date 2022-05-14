function out = MWIS(n,A,x)

%n = numnodes(G);
temp = de2bi(2^n - 1);
%A = adjacency(G);
max_sum = 0;
for i=1:2^n -1
    idx = logical(de2bi(i,size(temp,2)));
    a = A(idx,idx);
    if(a == zeros(size(a)) & sum(x(idx)) >= max_sum)
        max_sum = sum(x(idx));
        out = idx;
    end
end

end