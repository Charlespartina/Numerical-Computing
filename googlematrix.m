function [A,G,X,lambda] = googlematrix(M, alpha)
    %   Google Matrix Solver
    %   ------------------------------------------------
    %   INPUT:
    %       M:  A sparse matrix
    %       alpha:  A damping parameter
    %   ------------------------------------------------
    %   OUTPUT:
    %       A: The matrix A
    %       G:  The google matrix
    %       X: The eigenvectors of G
    %       lambda: The eigenvalues of G
    
    [m,n] = size(M);
    % Build matrix A
    A = zeros(n,n);
    col_sum = zeros(1,n);
    for i=1:n
        col_sum = col_sum + M(i,:);
    end
    for j=1:n
        if col_sum(1,j) == 0
            A(:,j) = 1/n;
        else
            A(:,j) = M(:,j)/col_sum(1,j);
        end
    end
    % Vector of all ones
    e = ones(n,1);
    % Compute matrix G using eig
    G = alpha*A + (1-alpha)/n*(e*e');
    [X,Lambda] = eig(G);
    lambda=diag(Lambda);
end