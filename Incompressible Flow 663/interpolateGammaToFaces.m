%% Now we find the gamma at the faces based on the Gfs
function [Gamma_e,Gamma_w,Gamma_n,Gamma_s] = interpolateGammaToFaces(Gamma,ge,gw,gn,gs,Nx,Ny)

Gamma_e = zeros(Nx,Ny);
Gamma_w = zeros(Nx,Ny);
Gamma_n = zeros(Nx,Ny);
Gamma_s = zeros(Nx,Ny);

for i = 1:Nx
    for j = 1:Ny
        if i < Nx
        Gamma_e(i,j)=1/((1-ge(i,j))/Gamma(i+1,j)+ge(i,j)/Gamma(i,j));
        end

        if i > 1
            Gamma_w(i,j)=1/((1-gw(i,j))/Gamma(i-1,j)+gw(i,j)/Gamma(i,j));
        end

        if j < Ny
            Gamma_n(i,j)=1/((1-gn(i,j))/Gamma(i,j+1)+gn(i,j)/Gamma(i,j));
      end

        if j > 1
            Gamma_s(i,j)=1/( (1-gs(i,j))/Gamma(i,j-1)+gs(i,j)/Gamma(i,j));
     end
    end
end
end