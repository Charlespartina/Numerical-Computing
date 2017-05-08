function [x,k,resi] = ConjugateGradient(A,b,x0,tol)
%   This function serves to solve linear equation Ax=b with CG method
%   Input:  A - The matrix A in the linear equation
%           b - The vector b in the linear equation
%           x0 - Initial guess for the solution
%           tol - Tolerance
%   Output: x - the solution
%           k - number of iteration
%           resi - residuals for each iteration
    x = x0;
    r = b-A*x;
    d= r'*r;
    resi = [];
    bd = b'*b;
    k = 1;
    p = r;
    while(d > (tol^2*bd) && k <= 200)
        resi(k) = sqrt(d);
        s = A*p;
        alpha=d/(p'*s);
        if(alpha <= 0)
            error('Alpha should be positive');
        end
        x=x+alpha*p;
        r=r-alpha*s;
        d1=r'*r;
        p=r+d1/d*p;
        d=d1;
        k=k+1;
    end
end