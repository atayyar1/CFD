function [aC, bC] = assembleTransientTerm(T_old, Vc,rho,cp,dt,aC,bC,Nx, Ny)
% adds the transient term to the system
% basically adds rho*cp*V/dt to diagonal and rho*cp*V/dt * T_old to RHS
for i = 1:Nx
    for j = 1:Ny
        transient_coeff = rho * cp * Vc(i,j) / dt;
        aC(i,j) = aC(i,j) +transient_coeff; bC(i,j)=bC(i,j)+transient_coeff * T_old(i,j);
    end
end
end