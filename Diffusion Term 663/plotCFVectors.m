function plotCFVectors(X,Y,XC,YC,CE,CW,CN,CS,Nx,Ny)

figure
hold on
axis equal
box on

%% Plot grid
[m,n] = size(X);

gridHandle = plot(nan,nan,'k','DisplayName','Grid lines');

for i = 1:m
    plot(X(i,:),Y(i,:),'k')
end

for j = 1:n
    plot(X(:,j),Y(:,j),'k')
end

%% Plot centroids
centroidHandle = plot(XC(:),YC(:),'bo','MarkerSize',6,'LineWidth',1.5, ...
    'DisplayName','Cell centroids');

%% Plot CF vectors
scale = 0.4;

ceHandle = plot(nan,nan,'r-','LineWidth',1.5,'DisplayName','CE: center to east neighbor');
cwHandle = plot(nan,nan,'g-','LineWidth',1.5,'DisplayName','CW: center to west neighbor');
cnHandle = plot(nan,nan,'b-','LineWidth',1.5,'DisplayName','CN: center to north neighbor');
csHandle = plot(nan,nan,'m-','LineWidth',1.5,'DisplayName','CS: center to south neighbor');

for i = 1:Nx
for j = 1:Ny

    xc = XC(i,j);
    yc = YC(i,j);

    if i < Nx
        quiver(xc,yc,scale*CE(i,j,1),scale*CE(i,j,2),'r')
    end

    if i > 1
        quiver(xc,yc,scale*CW(i,j,1),scale*CW(i,j,2),'g')
    end

    if j < Ny
        quiver(xc,yc,scale*CN(i,j,1),scale*CN(i,j,2),'b')
    end

    if j > 1
        quiver(xc,yc,scale*CS(i,j,1),scale*CS(i,j,2),'m')
    end

end
end

title('Center-to-Center Vectors (CF)')
xlabel('x')
ylabel('y')

legend([gridHandle,centroidHandle,ceHandle,cwHandle,cnHandle,csHandle], ...
    'Location','eastoutside')

hold off

end
