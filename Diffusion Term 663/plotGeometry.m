function plotGeometry(X,Y,XC,YC,Se,Sw,Sn,Ss,Nx,Ny, Xe, Ye,Xn,Yn,CE)

figure
hold on
axis equal
box on

%% ---- Plot grid lines ----
[m,n] = size(X);

gridHandle = plot(nan,nan,'k','DisplayName','Grid lines');

for i = 1:m
    plot(X(i,:),Y(i,:),'k')
end

for j = 1:n
    plot(X(:,j),Y(:,j),'k')
end

%% ---- Plot nodes ----
nodeHandle = plot(X(:),Y(:),'r.','MarkerSize',10,'DisplayName','Nodes');

%% ---- Plot centroids ----
centroidHandle = plot(XC(:),YC(:),'bo','MarkerSize',6,'LineWidth',1.5, ...
    'DisplayName','Cell centroids');

%% ---- Plot surface vectors ----
scale = 0.3;   % scale for visibility

SeHandle = plot(nan,nan,'r-','LineWidth',1.5,'DisplayName','S_e: east surface vector');
SwHandle = plot(nan,nan,'g-','LineWidth',1.5,'DisplayName','S_w: west surface vector');
SnHandle = plot(nan,nan,'b-','LineWidth',1.5,'DisplayName','S_n: north surface vector');
SsHandle = plot(nan,nan,'m-','LineWidth',1.5,'DisplayName','S_s: south surface vector');

for i = 1:Nx
    for j = 1:Ny

        xc = XC(i,j);
        yc = YC(i,j);

        quiver(xc,yc,scale*Se(i,j,1),scale*Se(i,j,2),'r')
        quiver(xc,yc,scale*Sw(i,j,1),scale*Sw(i,j,2),'g')
        quiver(xc,yc,scale*Sn(i,j,1),scale*Sn(i,j,2),'b')
        quiver(xc,yc,scale*Ss(i,j,1),scale*Ss(i,j,2),'m')

    end
end
scale = 0.4;

ceHandle = plot(nan,nan,'k-','LineWidth',1.5,'DisplayName','CE: center-to-east vector');

for i=1:Nx
    for j=1:Ny
        quiver(XC(i,j),YC(i,j),scale*CE(i,j,1),scale*CE(i,j,2),'k')
    end
end
xeHandle = plot(Xe(:),Ye(:),'gs','MarkerSize',6,'DisplayName','East face centers');
xnHandle = plot(Xn(:),Yn(:),'ms','MarkerSize',6,'DisplayName','North face centers');
title('Grid Geometry and Surface Vectors')
xlabel('x')
ylabel('y')
 
legend([gridHandle,nodeHandle,centroidHandle, ...
    SeHandle,SwHandle,SnHandle,SsHandle, ...
    ceHandle,xeHandle,xnHandle], ...
    'Location','eastoutside')

hold off

end
