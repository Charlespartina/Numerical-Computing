function [r,iterates,info] = findRoot(a,b,f,fderiv)
%function [r,iterates,info, iterates_diff, M] = findRoot(a,b,f,fderiv,fsecond)
    % Find a root of function f using hybrid of bisection method and
    % Newton's method.
    % --------------------------------------------------------------
    % INPUT:    a - lower bound of the root
    %           b - upper bound of the root
    %           f - the function
    %           fderiv - the derivative of f
    % OUTPUT:   r - the root
    %           iterates - all iterates of x_k
    %           info - 0 if the kth step if bisection and 1 if Newton
    
    iterates = zeros(0);
    % iterates_diff = zeros(0);
    info = [];
    % M=NaN;
    % Start with evaluating f(a) and f(b)
    fa = f(a);
    fb = f(b);
    iter = 1;
    if fa == 0
        r = a;
        return ;
    end
    if fb == 0
        r = b;
        return ;
    end
    if fa*fb > 0
        r = NaN;
        return ;
    end
    eps = 1e-13; % The machine epsilon
    % Loop
    x = a;
    r = a;
    fx = f(x);
    fpx = fderiv(x);
    while iter <= 200
        xnew = x-fx/fpx;
        if xnew >= a && xnew <= b
            info = [info ' N'];
        else
            xnew = (a+b)/2;
            info = [info ' B'];
        end
        iterates(iter) = xnew;
        fx = f(xnew);
        fpx = fderiv(xnew);
        if(fx*f(a) >= 0)
            a = xnew;
        else
            b = xnew;
        end
        if iter >= 2 && abs(iterates(iter)-iterates(iter-1))/max(1,abs(iterates(iter-1))) <= eps
            r = iterates(iter);
            break;
        end
        x = xnew;
        iter = iter + 1;
    end
%     % Convergence Ratio
%     iterates1 = [iterates 0];
%     iterates2 = [0 iterates];
%     iterates_diff = abs(iterates1-r)./abs(iterates2-r).^2;
%     [~,n] = size(iterates_diff);
%     iterates_diff = iterates_diff(2:n-1);
%     M = abs(fsecond(r))/(2*abs(fderiv(r)));
end