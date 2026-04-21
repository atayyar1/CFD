function [X,Y]=GridRect(m,n,Lx,Ly)
[X,Y] = meshgrid(linspace(0,Lx,m), linspace(0,Ly,n));
X = X'; 


Y =Y';
end