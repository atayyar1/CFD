function [X,Y] = Grid_test(m,n)

% Uniform square domain [0,1] × [0,1]

x = linspace(0,1,m);
y = linspace(0,1,n);

X = zeros(m,n);
Y = zeros(m,n);

for i = 1:m
    for j = 1:n
        X(i,j) = x(i);
        Y(i,j) = y(j);
    end
end

end