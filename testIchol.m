function [] = testIchol(A)
    alltol=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1];
    nz = zeros(7);
    dis = zeros(7);
    for i=1:7
        L=ichol(A,struct('type','ict','droptol',alltol(i)));
        nz(i)=nnz(L);
        dis(i)=norm(L*L'-A,1);
    end
    plot(log10(alltol), log10(nz), log10(alltol), log10(dis),'r-');
    ylabel('log10(Number of non zero entries) & log10(Norm of Discrepancy)');
    xlabel('log10(Tolerance)');
    legend('log(Number of non zero entries)', 'log(Norm of Discrepancy)');
end