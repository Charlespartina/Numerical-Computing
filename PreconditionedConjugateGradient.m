function [x,k,resi] = PreconditionedConjugateGradient(Avmult, b, x0, pi, tol)
    x = x0;
    r = b-Avmult(x);
    h = pi\r;
    d= r'*h;
    resi = [];
    bd = b'*(pi\b);
    k = 1;
    p = h;
    while(d > (tol^2*bd) && k <= 1000)
        resi(k) = sqrt(d);
        s = Avmult(p);
        alpha=d/(p'*s);
        if(alpha <= 0)
            error('Alpha should be positive');
        end
        x=x+alpha*p;
        r=r-alpha*s;
        h = pi\r;
        d1=r'*h;
        p=h+d1/d*p;
        d=d1;
        k=k+1;
    end
end