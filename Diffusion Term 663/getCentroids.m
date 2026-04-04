%% 1st Function: Centroid, Simple
function [XC, YC] = getCentroids(X, Y)
[m,n] = size(X);
XC = zeros(m-1,n-1);
YC = zeros(m-1,n-1);
for i = 1:m-1
    for j = 1:n-1
        XC(i,j) = 0.25*( X(i,j)   + X(i+1,j) + X(i,j+1) + X(i+1,j+1) );
        YC(i,j) = 0.25*( Y(i,j)   + Y(i+1,j) + Y(i,j+1) + Y(i+1,j+1) );

    end
end

end
