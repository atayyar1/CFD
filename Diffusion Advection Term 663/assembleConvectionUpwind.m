function [aE,aW,aN,aS,aC] = assembleConvectionUpwind(aE,aW,aN,aS,aC,Fe,Fw,Fn,Fs,Nx,Ny)
% Adds upwind convection contributions to diffusion coefficients, that were
% done in the old assignment
% Sign convention matches computeDiffusionCoefficients

for i = 1:Nx
    for j = 1:Ny
        % East face
        aE(i,j)=aE(i,j) - max(-Fe(i,j), 0);
        % West face
        aW(i,j)=aW(i,j)- max(-Fw(i,j),0);
        % North face
        aN(i,j) = aN(i,j) - max(-Fn(i,j), 0);
        
        % South face
        aS(i,j) = aS(i,j) - max(-Fs(i,j), 0);

      
        aC(i,j) = aC(i,j) + max(Fe(i,j),0) + max(Fw(i,j),0) + max(Fn(i,j),0) + max(Fs(i,j),0);
    end
end
end