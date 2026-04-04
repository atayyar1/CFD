function [X,Y] = Grid_test(m,n)
% Uniform square domain [0,1] x [0,1]
[X,Y] = meshgrid(linspace(0,1,m), linspace(0,1,n));
X = X'; Y = Y';
end