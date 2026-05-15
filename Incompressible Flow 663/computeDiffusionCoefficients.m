function [aE,aW,aN,aS,aC,bC] = computeDiffusionCoefficients(Gamma_e,Gamma_w,Gamma_n,Gamma_s,Ee,Ew,En,Es,CE,CW,CN,CS, ...
    S,Vc,bCorr,Nx,Ny)

aE = zeros(Nx,Ny);
aW = zeros(Nx,Ny);
aN = zeros(Nx,Ny);
aS = zeros(Nx,Ny);
aC = zeros(Nx,Ny);
bC = zeros(Nx,Ny);

for i = 1:Nx
    for j = 1:Ny
        % East
        if i < Nx
            aE(i,j) = -Gamma_e(i,j) * norm(squeeze(Ee(i,j,:))) / norm(squeeze(CE(i,j,:)));
        end
        %West
        if i > 1
            aW(i,j) = -Gamma_w(i,j) * norm(squeeze(Ew(i,j,:))) / norm(squeeze(CW(i,j,:)));
        end

        %North
        if j < Ny
            aN(i,j) = -Gamma_n(i,j) * norm(squeeze(En(i,j,:))) / norm(squeeze(CN(i,j,:)));
        end

        %% South
        if j > 1
            aS(i,j) = -Gamma_s(i,j) * norm(squeeze(Es(i,j,:))) / norm(squeeze(CS(i,j,:)));
        end

        %% Central coefficient
        aC(i,j) = -(aE(i,j) + aW(i,j) + aN(i,j) + aS(i,j));

        %% Right-hand side
        bC(i,j) = S(i,j)*Vc(i,j) + bCorr(i,j);

    end
end

end