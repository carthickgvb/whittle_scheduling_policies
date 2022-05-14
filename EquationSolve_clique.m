%Calculate Whittle Index for NON graphical Whittle Index based algorithm
function out = EquationSolve_clique(M,psi,C,l)

N = 500; %number of iterations
gamma = 0.1;

step = ceil(M/10);
x_iter = 0:step:M;

whittle_idx = zeros(size(x_iter,2),1);
X_axis = zeros(size(x_iter,2),1);
%tic;
idx=1;

%Solving the system of equation by converting to type AX = B
for x = x_iter
    lambda = 10*ones(N+1,1);
    
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
    end

    for y = x:M
        k = 0:M+min(y,psi)-y-1;
        A(y+2,y - min(y,psi) + k + 1) = mu(k,l);
        A(y+2,M+2) = -1;
        A(y+2,M+1) = -sum(A(y+2,:));
        A(y+2,y+1) = A(y+2,y+1) - 1; %-1 for LHS of our system of equations
        B(y+2,1) = -(C*y + f(min(y,psi)));
    end
    
    for n = 1:N
        y = 0:x-1;
        B(y+2,1) = -(C*y + lambda(n));
        
        X = linsolve(A,B);

        k = 0:M;
        temp_sum = mu(k,l)*(X(min(x - min(x,psi) + k,M) + 1) - X(min(x + k,M) + 1));
        lambda_old = lambda(n);
        lambda(n+1) = lambda(n) + gamma*(f(min(x,psi)) - lambda(n) + temp_sum);
        
        if (abs(lambda_old - lambda(n+1)) < 0.0001)
            break;
        end
    end
    
    whittle_idx(idx) = lambda(n+1);
    X_axis(idx) = x+1;
    idx=idx+1;    
end

    out = interp1(X_axis,whittle_idx,1:M+1,'linear','extrap');
end