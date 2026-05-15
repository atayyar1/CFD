function [aE, aW, aN, aS, aC, bC] = assemblePressureCorrection(Fe, Fw, Fn, Fs, aP_u, aP_v, Vc, rho, ...
    ge, gn, Nx, Ny, dx, dy, j_outlet)
aE = zeros(Nx,Ny);
aW = zeros(Nx,Ny);
aN = zeros(Nx,Ny);
aS = zeros(Nx,Ny);
aC = zeros(Nx,Ny);
bC = zeros(Nx,Ny);

 Du = Vc ./ aP_u;    
Dv = Vc ./ aP_v;    
 
for i = 1:Nx-1
    for j = 1:Ny
        gC   = ge(i,j);
        gF   = 1 - gC;
        D_e  = gC*Du(i,j) + gF*Du(i+1,j);   
        a_e  = -rho * D_e * dy / dx;     
        aE(i,j)   = aE(i,j)   + a_e;
        aW(i+1,j) = aW(i+1,j) + a_e;    
    end
end
 
for i = 1:Nx
    for j = 1:Ny-1
        gC   = gn(i,j);
        gF   = 1 - gC;
        D_n  = gC*Dv(i,j) + gF*Dv(i,j+1);
        a_n  = -rho * D_n * dx / dy;
        aN(i,j)   = aN(i,j)   + a_n;
        aS(i,j+1) = aS(i,j+1) + a_n;
    end
end
 
for j = 1:Ny
    if any(j_outlet == j)
        D_b      = Du(1,j);
        dist_Cb  = dx / 2;
        aC(1,j)  = aC(1,j) + rho * D_b * dy / dist_Cb;    end
end
aC = aC - (aE + aW + aN + aS); 
bC = -(Fe + Fw + Fn + Fs);
 if any(aC(:) <= 0)
    warning('assemblePressureCorrection: error')
end

end