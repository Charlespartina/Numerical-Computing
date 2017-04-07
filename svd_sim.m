function [] = svd_sim(num)
    clf;
    % Generate points on unit circle
    % Each column of v is a point on the unit circle
    v = zeros(2,num);
    ls = linspace(0,2*pi, num);
    for i=1:length(ls)
        [v(1,i), v(2,i)] = generate_v(ls(i));
    end
    
    % Create random matrix
    A = randn(2);
    % Transform v into u
    u = A*v;
    % Compute SVD
    [U,D,V] = svd(A)
    
    % Plot unit circle
    subplot(1,2,1)
    plot(v(1,:),v(2,:));
    hold on
    % Plot v1 and v2 by retreving columns of V
    plot([0,V(1,1)], [0,V(2,1)]);
    plot([0,V(1,2)], [0,V(2,2)]);
    % Adjust the plot
    title('v')
    grid on;
    daspect([1 1 1]);
    
    % Plot ellipse
    subplot(1,2,2)
    plot(u(1,:),u(2,:));
    hold on
    % Plot u1 and u2 by retreving columns of U and D
    plot([0,U(1,1)*D(1,1)], [0,U(2,1)*D(1,1)]);
    plot([0,U(1,2)*D(2,2)], [0,U(2,2)*D(2,2)]);
    % Adjust the plot
    title('u')
    grid on;
    daspect([1 1 1]);
end

function [x,y] = generate_v(alpha)
    x = cos(alpha);
    y = sin(alpha);
end