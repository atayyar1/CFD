%% for the interior faces! s
function bCorr = computeNonOrthogonalCorrection(Gamma_e,Gamma_w,Gamma_n,Gamma_s,dTdx,dTdy,Tev,Twv,Tnv,Tsv,ge,gw,gn,gs,Nx,Ny)

bCorr = zeros(Nx,Ny);

for i = 1:Nx
    for j = 1:Ny

        corr = 0;

        %% East fac
        if i < Nx
            dTdx_f = ge(i,j)*dTdx(i,j) + (1-ge(i,j))*dTdx(i+1,j);
            dTdy_f = ge(i,j)*dTdy(i,j) + (1-ge(i,j))*dTdy(i+1,j);
            gradTf = dTdx_f*Tev(i,j,1) + dTdy_f*Tev(i,j,2);
            corr = corr + Gamma_e(i,j)*gradTf;
        end

        %% West face
        if i > 1
            dTdx_f = gw(i,j)*dTdx(i,j) + (1-gw(i,j))*dTdx(i-1,j);
            dTdy_f = gw(i,j)*dTdy(i,j) + (1-gw(i,j))*dTdy(i-1,j);
            gradTf = dTdx_f*Twv(i,j,1) + dTdy_f*Twv(i,j,2);
            corr = corr + Gamma_w(i,j)*gradTf;
        end

        %% North face
        if j < Ny
            dTdx_f = gn(i,j)*dTdx(i,j) + (1-gn(i,j))*dTdx(i,j+1);
            dTdy_f = gn(i,j)*dTdy(i,j) + (1-gn(i,j))*dTdy(i,j+1);
            gradTf = dTdx_f*Tnv(i,j,1) + dTdy_f*Tnv(i,j,2);
            corr = corr + Gamma_n(i,j)*gradTf;
        end

        %% South face
        if j > 1
            dTdx_f = gs(i,j)*dTdx(i,j) + (1-gs(i,j))*dTdx(i,j-1);
            dTdy_f = gs(i,j)*dTdy(i,j) + (1-gs(i,j))*dTdy(i,j-1);
            gradTf = dTdx_f*Tsv(i,j,1) + dTdy_f*Tsv(i,j,2);
            corr = corr + Gamma_s(i,j)*gradTf;
        end

        bCorr(i,j) = corr;

    end
end

end
