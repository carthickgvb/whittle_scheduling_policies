function out = MWIS_neelay(n,A,x,psi)

%n = numnodes(G);
temp = de2bi(2^n - 1);
%A = adjacency(G);
max_sum = Inf*(-1);
v1=200;
    conc=horzcat(x,psi);
    Z=min(conc,[],2);
    first=x.*Z;
    V=v1*ones(n,1);
    Fun=zeros(n,1);
    for k=1:n
        Fun(k)=f(Z(k));
    end
for i=1:2^n -1
    idx = logical(de2bi(i,size(temp,2)));
    a = A(idx,idx);
    %W=x.*transpose(idx);

    
   % Fun=arrayfun(f,Z);
    Rig=V.*Fun;
    result=first-Rig;
    if(a == zeros(size(a)) & sum(result(idx)) >= max_sum)
        max_sum = sum(result(idx));
        out = idx;
    end
end