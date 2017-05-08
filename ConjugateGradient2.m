function [x,k,resi] = ConjugateGradient2(Avmult,b,x0,tol)
%   This function serves to solve linear equation Ax=b with CG method
%   Input:  Avmult - A function that computes A*v
%           b - The vector b in the linear equation
%           x0 - Initial guess for the solution
%           tol - Tolerance
%   Output: x - the solution
%           k - number of iteration
%           resi - residuals for each iteration
    x = x0;
    r = b-Avmult(x);
    d= r'*r;
    resi = [];
    bd = b'*b;
    k = 1;
    p = r;
    while(d > (tol^2*bd) && k <= 1000)
        resi(k) = sqrt(d);
        s = Avmult(p);
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