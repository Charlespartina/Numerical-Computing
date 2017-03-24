function [A,x,p] = polyInterpOrApprox(t,b,deg,nPlot)
%   Polynomial Interpolation or Approximation
%   ------------------------------------------------
%   INPUT:
%       t:  data vector which is ordered and distinct
%       b:  data vector with any real entries
%       deg:    degree of polynomial for interpolation
%       nPlot:  the number of equally spaced grid of points
%   ------------------------------------------------
%   OUTPUT:
%       A:  Vandermonde matrix for interpolating/approximating the data
%       x:  coefficient for polynomial
%       p:  polynomial values

m = 0; % The number of data input
% Check the dimensions of t and b
if(length(t) ~= length(b))
    error('Data vector t and data vector b must have same length');
else
    m = length(t);
end
if(m < deg+1)
    error('Not enough data records for determining the polynomial');
end
powers = 0:deg;  % Degrees of polynomial
A = zeros(m,deg+1); % A zero matrix
for j=1:deg+1
    A(:,j) = t.^powers(j);
end
x = A \ b; % Determine the coefficient for the polynomial

tt = linspace(t(1), t(m), nPlot); % A sequence of equally spaced points 
p = zeros(1,nPlot);
for i = 1:deg+1
    p = p + x(i)*tt.^(i-1);
end
plot(tt,p);
hold on;
plot(t',b','ro','LineWidth',2);
xlabel('t');
ylabel('v');
hold off;
end