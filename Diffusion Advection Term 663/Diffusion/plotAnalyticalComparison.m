function plotAnalyticalComparison(Tsol,T_exact,errorField,XC,YC,m,n)
%PLOTANALYTICALCOMPARISON Compare the analytical and numerical benchmark solutions.

figure
tiledlayout(1,3,'TileSpacing','compact','Padding','compact')

nexttile
contourf(XC,YC,T_exact,20,'LineColor','none')
colorbar
xlabel('X')
ylabel('Y')
title(sprintf('Analytical (%d x %d)',m,n))

nexttile
contourf(XC,YC,Tsol,20,'LineColor','none')
colorbar
xlabel('X')
ylabel('Y')
title(sprintf('Numerical (%d x %d)',m,n))

nexttile
contourf(XC,YC,errorField,20,'LineColor','none')
colorbar
xlabel('X')
ylabel('Y')
title('|T_{num} - T_{exact}|')

sgtitle('Uniform-Grid Benchmark: Theoretical vs Numerical Solution')

end
