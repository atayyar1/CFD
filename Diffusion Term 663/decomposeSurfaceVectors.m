function [Ee,Ew,En,Es,Tev,Twv,Tnv,Tsv] = decomposeSurfaceVectors(Se,Sw,Sn,Ss,CE,CW,CN,CS,Nx,Ny,method)

Ee  = zeros(Nx,Ny,2);
Ew  = zeros(Nx,Ny,2);
En  = zeros(Nx,Ny,2);
Es  = zeros(Nx,Ny,2);

Tev = zeros(Nx,Ny,2);
Twv = zeros(Nx,Ny,2);
Tnv = zeros(Nx,Ny,2);
Tsv = zeros(Nx,Ny,2);

for i = 1:Nx
    for j = 1:Ny

        %% East
        if i < Nx
            sf = [Se(i,j,1) Se(i,j,2)];
            cf = [CE(i,j,1) CE(i,j,2)];
            [Ef,Tf] = oneFaceDecomposition(sf,cf,method);
            Ee(i,j,:)  = Ef;
            Tev(i,j,:) = Tf;
        end

        %% West
        if i > 1
            sf = [Sw(i,j,1) Sw(i,j,2)];
            cf = [CW(i,j,1) CW(i,j,2)];
            [Ef,Tf] = oneFaceDecomposition(sf,cf,method);
            Ew(i,j,:)  = Ef;
            Twv(i,j,:) = Tf;
        end

        %% North
        if j < Ny
            sf = [Sn(i,j,1) Sn(i,j,2)];
            cf = [CN(i,j,1) CN(i,j,2)];
            [Ef,Tf] = oneFaceDecomposition(sf,cf,method);
            En(i,j,:)  = Ef;
            Tnv(i,j,:) = Tf;
        end

        %% South
        if j > 1
            sf = [Ss(i,j,1) Ss(i,j,2)];
            cf = [CS(i,j,1) CS(i,j,2)];
            [Ef,Tf] = oneFaceDecomposition(sf,cf,method);
            Es(i,j,:)  = Ef;
            Tsv(i,j,:) = Tf;
        end

    end
end

end

function [Ef,Tf] = oneFaceDecomposition(sf,cf,method)

ef = cf / norm(cf);              % unit vector along CF
Sfmag = norm(sf);
n = sf / Sfmag;                  % unit normal direction of the face
cosTheta = dot(n,ef);

if strcmp(method,'minimum')
    Ef = dot(sf,ef) * ef;

elseif strcmp(method,'normal')
    Ef = Sfmag * ef;

elseif strcmp(method,'overrelaxed')
    Ef = (Sfmag / cosTheta) * ef;

else
    error('Unknown method. Use minimum, normal, or overrelaxed.');
end

Tf = sf - Ef;

end