%% Construct Matrix for general problem
function [A,x] = TrussGeneral(k, load, op)
%   Generalized Truss Problem Solver
%   --------------------------------------------------------------
%   INPUT:
%       k:  Number of sections
%       load:  The loads under intersections
%   --------------------------------------------------------------
%   OUTPUT:
%       x:  The internal force for each bars

    % Check if the load vector has correct size
    if length(load) ~= 2*k+1
        error('The size of load vector is inconsistent with k');
    end
    % Dimension of matrix and extended load vector.
    n = 2+3+k*8;
    A = zeros(n, n);
    b=zeros(n,1); % Extended load vector
    alpha = cos(0.25*pi);
    for i=1:k
        % Starting index of this section
        % We let each loop be responsible for 8 equations from 4 joints
        offsetC = 2+(i-1)*8;
        offsetR = (i-1)*8;
        % Two equations for upper-left joint
        % If the section is the left-most section, then we have to handle
        % the equation specially
        if i == 1
            A(1,1) = alpha;
            A(1,4) = -1;
            A(1,5) = -alpha;
            A(2,3) = 1;
            A(2,1) = alpha;
            A(2,5) = alpha;
        else
            A(offsetR+1, offsetC+1-3) = 1;
            A(offsetR+1, offsetC+1+1) = -1;
            A(offsetR+1, offsetC+1-2) = alpha;
            A(offsetR+1, offsetC+1+2) = -alpha;
            A(offsetR+2, offsetC+1) = 1;
            A(offsetR+2, offsetC+1+2) = alpha;
            A(offsetR+2, offsetC+1-2) = alpha;
        end
        % Two equations for lower-left joint
        A(offsetR+3, offsetC+1-1) = 1;
        A(offsetR+3, offsetC+1+3) = -1;
        A(offsetR+4, offsetC+1) = 1;
        % Two equations for upper-middle joint
        A(offsetR+5, offsetC+1+1) = 1;
        A(offsetR+5, offsetC+1+5) = -1;
        A(offsetR+6, offsetC+1+4) = 1;
        % Two equations for lower-middle joint
        A(offsetR+7, offsetC+1+3) = 1;
        A(offsetR+7, offsetC+1+2) = alpha;
        A(offsetR+7, offsetC+1+6) = -alpha;
        A(offsetR+7, offsetC+1+7) = -1;
        A(offsetR+8, offsetC+1+4) = 1;
        A(offsetR+8, offsetC+1+2) = alpha;
        A(offsetR+8, offsetC+1+6) = alpha;
    end
    % Build equations for the last two joints and the right moving part
    A(n-4,n-5) = 1;
    A(n-4,n-4) = alpha;
    A(n-4,n-1) = -alpha;
    A(n-3,n-4) = alpha;
    A(n-3,n-2) = 1;
    A(n-3,n-1) = alpha;
    A(n-2,n-3) = 1;
    A(n-2,n) = -1;
    A(n-1,n-2) = 1;
    A(n,n-1) = alpha;
    A(n,n) = 1;
    % Build vector b
    for k=1:length(load)
        b(k*4,1) = load(k);
    end
    A = sparse(A);
    if op == 1
        x = A\b;
    elseif op ==2
        x = inv(A)*b;
    else
        [l,u] = lu(A);
        y = l\b;
        x = u\y;
    end
end