function [A,kA,x1,x2,x3,x4] = approxpoly(t,b,d,cut)
%   Polynomial Interpolation or Approximation
%   Solve coefficient with 4 different methods.
%   ------------------------------------------------
%   INPUT:
%       t:  data vector which is ordered and distinct
%       b:  data vector with any real entries
%       d:    degree of polynomial for interpolation
%       cut:  cutoff tolerance for SVD method
%   ------------------------------------------------
%   OUTPUT:
%       A:  Vandermonde matrix for interpolating/approximating the data
%       kA: condition number of A
%       x1: coefficient for polynomial solved by using backslash
%       x2: coefficient for polynomial solved by normal equation
%       x3: coefficient for polynomial solved by QR decomposition
%       x4: coefficient for polynomial solved by SVD
    
    if(length(t) ~= length(b))
        error('Data vector t and data vector b must have same length');
    else  
        m = length(t);
    end

    if(m >= d+1)
        n = d+1;
        powers = 0:d;  % Degrees of polynomial
        A = zeros(m,d+1); % A zero matrix
        for j=1:d+1
            A(:,j) = t.^powers(j);
        end
        % Different Method for computing 
        % 1. Using Backslash
        x1 = A\b;
         
        % 2. Normal Equations
        B = A'*A;
        y = A'*b;
        [L,p] = chol(B, 'lower');
        if p ~= 0
            x2 = NaN(n,1);
        else
            z = forsub(L,y);
            x2 = backsub(L',z);
        end
        
        % 3. QR Decomposition
        [Q,R] = qr(A,0); % Economy QR
        c = Q'*b;
        x3 = backsub(R,c);
        
        % 4. SVD
        [U,S,V] = svd(A,0); % Economy SVD
        z = U'*b;
        y = zeros(n,1);
        for i=1:n
            if S(i,i) > cut
                y(i) = z(i)/S(i,i);
            else
                y(i) = 0;
            end
        end
        x4=V*y;
        kA = S(1,1)/S(n,n);
    else
        error('m must be greater or equal to d+1');
    end
    
    
    
end

%% Function for Forward Substitution
function y = forsub(A,b)
    n = length(b);
    % Make unit lower triangular matrix
    for i = 1:n
        dig = A(i,i);
        A(i,:) = A(i,:)./dig;
        b(i) = b(i)./dig;
    end
    y = zeros(n,1);
    y(1) = b(1);
    for k = 2:n
        y(k) = b(k) - A(k, 1:k-1) * y(1:k-1);
    end
end

%% Function for Backward Substitution
function x = backsub(A,b)
    n = length(b);
    x = b;
    x(n) = b(n)/A(n,n);
    for k = n-1:-1:1
        x(k) = (b(k)-A(k,k+1:n)*x(k+1:n))/A(k,k);
    end
end