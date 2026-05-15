function [aE, aW, aN, aS, aC, bC, aP] = assembleMomentumCoeffs( ...
    comp, phi, rho, mu, p, ...
    Fe, Fw, Fn, Fs, ...
    Ee, Ew, En, Es, CE, CW, CN, CS, ...
    Se, Sw, Sn, Ss, Vc, Nx, Ny, alpha)

mu_e = mu *ones(Nx, Ny);
mu_w = mu *ones(Nx, Ny);
mu_n = mu *ones(Nx, Ny);
mu_s = mu *ones(Nx, Ny);

S_zero = zeros(Nx, Ny);
bCorr_zero=zeros(Nx, Ny);

[aE, aW, aN, aS, aC, bC] = computeDiffusionCoefficients(mu_e, mu_w, mu_n, mu_s,Ee, Ew, En, Es, CE, CW, CN, CS,S_zero, Vc, bCorr_zero, Nx, Ny); 
[aE, aW, aN, aS, aC] = assembleConvectionUpwind(aE, aW, aN, aS, aC, Fe, Fw, Fn, Fs, Nx, Ny);
aP = aC;  
for i = 1:Nx
    for j = 1:Ny
    if i < Nx; p_e = 0.5*(p(i,j) + p(i+1,j));  else;  p_e = p(Nx,j);  end
        if i > 1;   p_w = 0.5*(p(i-1,j) + p(i,j));  else;  p_w = p(1,j);   end
       if j < Ny;  p_n = 0.5*(p(i,j) +   p(i,j+1));  else;  p_n = p(i,Ny);  end
        if j > 1;   p_s = 0.5*(p(i,j-1) + p(i,j));  else;  p_s = p(i,1);   end
      switch lower(comp)
           case 'u'
                area_x = norm(squeeze(Se(i,j,:)));
                bC(i,j) = bC(i,j) - (p_e - p_w) * area_x;
            case 'v'
               area_y =   norm(squeeze(Sn(i,j,:)));
                bC(i,j) = bC(i,j) - (p_n - p_s) * area_y;
            otherwise
                error('assembleMomentumCoeffs: error', comp)
      end
    end
end

if any(aC(:) <= 0)
   warning('assembleMomentumCoeffs: non-positive diagonal', comp)
end
   end