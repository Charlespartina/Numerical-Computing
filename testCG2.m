function [] = testCG2(N, r, tol)
    G = numgrid('S',N+2);
    D = sprand(speye(N^2));
    K = delsq(G);
    B = sprand(N^2,r,1);
    Av = @(v)(D*v+K\(B*(B'*(K\v))));
    b = ones(N^2,1);
    x0 = zeros(N^2,1);
    [x,k,resi] = ConjugateGradient2(Av,b,x0,tol);
    x = reshape(x,N,N);
    figure()
    mesh(x);
    figure()
    semilogy(linspace(1,k-1,k-1), resi), title('Residual'), xlabel('k'), ylabel('log of residual')
end