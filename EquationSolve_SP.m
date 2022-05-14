function out = EquationSolve_SP(S,w_idx,G,M_bar,psi_bar,C,l_bar,type)
% psi = realmax;

n_nodes = length(S);
N = 100; %number of iterations
gamma = 0.01;
%type = 1;

% Solving the system of equation by converting to type AX = B
M_max = max(M_bar);
A_all = zeros(M_max+2,M_max+2,n_nodes);
B_all = zeros(M_max+2,n_nodes);

lambda_old = zeros(n_nodes,N+1);                                                        
lambda = zeros(n_nodes,1);   

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
        B(y+2,1) = -(C*y + sum(lambda(nbr)) + lambda(i));      %%
    end

    for y = x:M
        k = 0:M+min(y,psi)-y-1;
        A(y+2,y - min(y,psi) + k + 1) = mu(k,l);
        A(y+2,M+2) = -1;
        A(y+2,M+1) = -sum(A(y+2,:));
        A(y+2,y+1) = A(y+2,y+1)   - 1; %-1 for LHS of our system of equations
        B(y+2,1) = -(C*y + f(min(y,psi)));
    end
    
    A_all(1:M+2,1:M+2,i) = A;
    B_all(1:M+2,i) = B;
end

% iterating over lambdas
for n = 1:N
    for i = 1:n_nodes
        x = S(i);                                                          %%
        M = M_bar(i);                                                      %%
        psi = psi_bar(i);                                                  %%
        l = l_bar(i);                                                      %%

        X = linsolve(A_all(1:M+2,1:M+2,i),B_all(1:M+2,i));
        
        k = 0:M;
        temp_sum = mu(k,l)*(X(min(x - min(x,psi) + k,M) + 1) - X(min(x + k,M) + 1));
        if type==1
            lambda(i) = lambda_old(i,n) + gamma*(f(min(x,psi)) - lambda_old(i,n) + temp_sum);
        elseif type==2     
            lambda(i) = lambda_old(i,n) + gamma*(f(min(x,psi)) - sum(lambda_old(nbr,n)) - lambda_old(i,n) + temp_sum);
        end
    end
    
    for i=1:n_nodes
        x = S(i);
        nbr = neighbors(G,i);
        
        y = 0:x-1;
        B_all(y+2,i) = -(C*y + sum(lambda(nbr)) + lambda(i));
    end
    
    lambda_old(:,n+1) = lambda;
end
% for i = 1:n_nodes
%      plot(lambda_old(i,:));
%      hold on;
% end
% plot(lambda_old(1,:));
% hold on;
out = lambda;
end
