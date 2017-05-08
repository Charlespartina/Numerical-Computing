function [] = testPCG(N, r, tol)
    D = sprand(speye(N^2));
    G = numgrid('S',N+2);
    K = delsq(G);
    B = sprand(N^2,r,1);
    Av = @(v)(D*v+(K\(B*(B'*(K\v)))));
    b = ones(N^2,1);
    x0 = zeros(N^2,1);
    [x,k,resi] = PreconditionedConjugateGradient(Av,b,x0,D,tol);
    x = reshape(x,N,N);
    figure()
    mesh(x);
    figure()
    semilogy(linspace(1,k-1,k-1), resi), title('Residual'), xlabel('k'), ylabel('log of residual')
end