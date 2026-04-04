%% Taken from the Adam's Notes
function [X,Y] = Grid(m,n)

xi  = linspace(0,1,m);
eta = linspace(0,1,n);

X = zeros(m,n);
Y = zeros(m,n);

for i = 1:m
    Xi = xi(i);

    for j = 1:n
        Eta = eta(j);

        XY = (1-Eta)*Xb(Xi) + Eta*Xt(Xi) + (1-Xi)*Xl(Eta) + Xi*Xr(Eta) - (Xi*Eta*Xt(1) + Xi*(1-Eta)*Xb(1) + ...
              Eta*(1-Xi)*Xt(0) + (1-Xi)*(1-Eta)*Xb(0));

        X(i,j) = XY(1);
        Y(i,j) = XY(2);

    end
end

end

function p = Xb(s)

x = 0 + s*(2-0);
y = 0 + s*(0-0);

p = [x; y];

end


% Top:B-C
function p = Xt(s)

x = -1 + s*(3+1);
y = 3  + s*(5-3);

p = [x; y];

end


% Left: AB
function p = Xl(s)

x = 0  + s*(-1-0);
y = 0  + s*(3-0);

p = [x; y];

end


% Right: DC
function p = Xr(s)

x = 2 + s*(3-2);
y = 0 + s*(5-0);

p = [x; y];

end