%% Construct Matrix for original problem
function [A,x] = Truss()
%   Truss Problem Solver
%   --------------------------------------------------------------
%   INPUT:
%   --------------------------------------------------------------
%   OUTPUT:
%       x:  The internal force for each bars

    alpha = cos(0.25*pi);
    A=[
        alpha   0   0   -1  -alpha  0   0   0   0   0   0   0   0;
        alpha   0   1   0   alpha   0   0   0   0   0   0   0   0;
        0   1   0   0   0   -1  0   0   0   0   0   0   0;
        0   0   1   0   0   0   0   0   0   0   0   0   0;
        0   0   0   1   0   0   0   -1  0   0   0   0   0;
        0   0   0   0   0   0   1   0   0   0   0   0   0;
        0   0   0   0   alpha   1   0   0   -alpha  -1  0   0   0;
        0   0   0   0   alpha   0   1   0   alpha   0   0   0   0;
        0   0   0   0   0   0   0   1   alpha   0   0   -alpha  0;
        0   0   0   0   0   0   0   0   alpha   0   1   alpha   0;
        0   0   0   0   0   0   0   0   0   1   0   0   -1;
        0   0   0   0   0   0   0   0   0   0   1   0   0;
        0   0   0   0   0   0   0   0   0   0   0   alpha   1;
    ];
    b=[0;0;0;-20;0;0;0;0;0;0;0;20;0];
    x=A\b;
    draw(x);
end
%% Plot
function draw(x)
    % X axis of two ends of each bar
    pX=[
        [0,1];
        [0,1];
        [1,1];
        [1,2];
        [1,2];
        [1,2];
        [2,2];
        [2,3];
        [2,3];
        [2,3];
        [3,3];
        [3,4];
        [3,4];
    ];
    % Y axis of two ends of each bar
    pY=[
        [0,1];
        [0,0];
        [1,0];
        [1,1];
        [1,0];
        [0,0];
        [1,0];
        [1,1];
        [0,1];
        [0,0];
        [1,0];
        [1,0];
        [0,0];
    ];
    hold on;
    x2 = x/2; % x2 is a coly of x that has a smaller magnitude
    for i=1:13
        if x2(i) > 0+1e-9
            plot(pX(i,:),pY(i,:),'r-','LineWidth',x2(i));
        elseif x2(i) < 0-1e-9
            plot(pX(i,:),pY(i,:),'b-','LineWidth',-x2(i));
        else
            plot(pX(i,:),pY(i,:),'g--','LineWidth',1);
        end
    end
end