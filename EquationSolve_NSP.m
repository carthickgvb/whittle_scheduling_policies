function out = EquationSolve_NSP(S,w_idx,G,M_bar,psi_bar,C,l_bar,type)

n_nodes = length(S);
gamma = 0.001;
%type = 2;

% Solving the system of equation by converting to type AX = B
lambda = zeros(n_nodes,1);
%lambda = w_idx;

for i = 1:n_nodes
    x = S(i);                                                          %%
    M = M_bar(i);                                                      %%
    psi = psi_bar(i);                                                  %%
    l = l_bar(i);                                                      %%
    nbr = neighbors(G,i);                                              %%

    A = zeros(M+2,M+2);
    B = zeros(M+2,1);
    %matrix V(0) V(1) ... V(M) beta

    A(1,1) = -1; %V(0) = 0

    for y = 0:(x-1)
        k = 0:M-y-1;
        A(y+2,y+k+1) = mu(k,l);
        A(y+2,M+2) = -1;
        A(y+2,M+1) = - sum(A(y+2,:));
        A(y+2,y+1) = A(y+2,y+1) - 1; %-1 for LHS of our system of equations
        B(y+2,1) = -(C*y + sum(w_idx(nbr)) + w_idx(i));
        %B(y+2,1) = -(C*y + sum(lambda(nbr)) + lambda(i));
    end

    for y = x:M
        k = 0:M+min(y,psi)-y-1;
        A(y+2,y - min(y,psi) + k + 1) = mu(k,l);
        A(y+2,M+2) = -1;
        A(y+2,M+1) = -sum(A(y+2,:));
        A(y+2,y+1) = A(y+2,y+1)   - 1; %-1 for LHS of our system of equations
        B(y+2,1) = -(C*y + f(min(y,psi)));
    end
    X = linsolve(A,B);
    k = 0:M;
    temp_sum = mu(k,l)*(X(min(x - min(x,psi) + k,M) + 1) - X(min(x + k,M) + 1));
    if type==1
        lambda(i) = w_idx(i) + gamma*(f(min(x,psi)) - w_idx(i) + temp_sum);
        %lambda(i) = lambda(i) + gamma*(f(min(x,psi)) - lambda(i) + temp_sum);
    elseif type==2
        lambda(i) = w_idx(i) + gamma*(f(min(x,psi)) - sum(w_idx(nbr)) - w_idx(i) + temp_sum);
        %lambda(i) = lambda(i) + gamma*(f(min(x,psi)) - sum(lambda(nbr)) - lambda(i) + temp_sum);
    end
end

out = lambda;
end
