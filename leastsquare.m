function [] = leastsquare(t,b,d,nPlot)
    % Solve the least square problem and Plot
    [A,kA,x1,x2,x3,x4] = approxpoly(t',b,d,0);
    plotting(t',x4,b,d,nPlot);
    
    % Check residual for x4 and relative error for x1->x3 when d increases
    for dd=1:d
        % Record the residual for all methods
        [A,kA,x1,x2,x3,x4] = approxpoly(t',b,dd,0);
        resi1(dd) = norm(b-A*x1); % Compute the residual for Backslash
        resi2(dd) = norm(b-A*x2); % Compute the residual for Normal Equation
        resi3(dd) = norm(b-A*x3); % Compute the residual for QR
        resi4(dd) = norm(b-A*x4); % Compute the residual for SVD
        % Record the relative error for the first three methods
        eta1(dd) = norm(x1-x4)./norm(x4);
        eta2(dd) = norm(x2-x4)./norm(x4);
        eta3(dd) = norm(x3-x4)./norm(x4);
        condA(dd) = kA;
        % disp(condA');
    end
    deg = linspace(1,d,d);
    
    % Plot residual for x4 with different degree
    figure();
    plot(deg,resi4,'LineWidth',2);
    title('Residuals for SVD method');
    xlabel('Degree');
    ylabel('Residual');
    legend('SVD');
    
    % Plot residuals for all methods with different degree
    figure();
    plot(deg,resi1,'LineWidth',2);
    hold on;
    plot(deg,resi2,'LineWidth',2);
    plot(deg,resi3,'LineWidth',2);
    plot(deg,resi4,'LineWidth',2);
    hold off;
    title('Residuals for different methods');
    xlabel('Degree');
    ylabel('Residual');
    legend('Backslash','Normal Equation','QR Decomposition','SVD');
    
    % Plot relative errors
    figure();
    semilogy(deg,eta1);
    hold on;
    semilogy(deg,eta2);
    semilogy(deg,eta3);
    hold off;
    title('Relative Error with SVD solution');
    xlabel('Degree');
    ylabel('Relative Error (logarithmic scaled)');
    legend('Backslash','Normal Equation','QR Decomposition');
    
    % Plot condition numbers
    condB = condA.^2;
    % disp(condB');
    figure();
    semilogy(deg, condA);
    hold on;
    semilogy(deg, condB);
    hold off;
    title('Condition Numbers');
    xlabel('Degree');
    ylabel('Condition Number (logarithmic scaled)');
    legend('k(A)','k(B)');
    
    % Plot polynomial for both x2 and x4
    plotting2(t',x2,x4,b,d,nPlot);
    
    % Plot residuals and norm of x4 with different cutoff tolerance
    cuts = linspace(1,10,10);
    for cut=1:10
        [A,kA,x1,x2,x3,x4] = approxpoly(t',b,d,cut);
        resi_cut4(cut) = norm(b-A*x4); % Compute the residual for SVD
        norm4(cut) = norm(x4);
    end
    figure();
    plot(cuts,norm4,'LineWidth',2);
    hold on;
    plot(cuts,resi_cut4,'LineWidth',2);
    hold off;
    title('Residuals of SVD method and Norm of x4');
    xlabel('Cutoff tolerance');
    ylabel('Value');
    legend('Norm of x4','Residual');
end

%% Function for plotting the polynomial
function [] = plotting(t, x, b, d, nPlot)
    m = size(t);
    m = m(1,1);
    tt = linspace(t(1), t(m), nPlot); % A sequence of equally spaced points 
    p = zeros(1,nPlot);
    for i = 1:d+1
        p = p + x(i)*tt.^(i-1);
    end
    figure();
    plot(tt,p,'b-');
    hold on;
    plot(t',b','ro','LineWidth',2);
    xlabel('t');
    ylabel('v');
    title('Polynomial Interpolation/Approximation Result');
    legend('Polynomial','Input data')
end

%% Function for plotting the polynomial with both x2 and x4
function [] = plotting2(t, x2, x4, b, d, nPlot)
    m = size(t);
    m = m(1,1);
    tt = linspace(t(1), t(m), nPlot); % A sequence of equally spaced points 
    p2 = zeros(1,nPlot);
    p4 = zeros(1,nPlot);
    for i = 1:d+1
        p2 = p2 + x2(i)*tt.^(i-1);
        p4 = p4 + x4(i)*tt.^(i-1);
    end
    figure();
    plot(tt,p2,'b-','LineWidth',2);
    hold on;
    plot(tt,p4,'g-','LineWidth',2);
    plot(t',b','ro','LineWidth',2);
    xlabel('t');
    ylabel('v');
    title('Polynomial Interpolation/Approximation Result for both x2 and x4');
    legend('Polynomial for x2','Polynomial for x4','Input data')
end
