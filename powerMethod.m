%% Power Method
% function [v] = powerMethod(A, alpha, url)
%     n = size(A);
%     n = n(1);
%     e = ones(n,1);
%     v0 = (1/n)*e;
%     v1 = alpha*A*v0+(1-alpha)*(e'*v0)*e/n;
%     num = 100; % Number of iteration
%     mu = zeros(1,num);
%     for iter=1:num
%         v0 = v1;
%         % v1 = Gv0
%         v1 = alpha*A*v0+(1-alpha)*(e'*v0)*e/n;
%         mu(1,iter) = norm(v1-v0,1);
%     end
%     % Plot rate of convergence
%     [xsort, index] = sort(v0, 'descend');
%     
%     % Display the rank
%     for j=1:10
%         fprintf('index = %d ',index(j));
%         disp(url(index(j)));
%     end
%     
%     figure();
%     plot(linspace(1,num,num), mu);
%     title('Convergence Rate of Residuals');
%     xlabel('Num of iterations');
%     ylabel('Residual');
%     v = v0;
% end

 %% Blocked Power Method
function [V,T] = powerMethod(A,m,alpha,iter)
    n = size(A);
    n = n(1);
    e = ones(n,1);
    % Initial Orthogonal Matrix
    V = eye(n,m);
    V1 = V;
    % matrix-vector product
    for j=1:m
        V1(:,j) = alpha*A*V(:,j)+(1-alpha)*(e'*V(:,j))*e/n;
    end
    for iter=1:iter
        [V,R] = qr(V1,0);
        for j=1:m
            V1(:,j) = alpha*A*V(:,j)+(1-alpha)*(e'*V(:,j))*e/n;
        end
        T = V'*V1;
    end
end