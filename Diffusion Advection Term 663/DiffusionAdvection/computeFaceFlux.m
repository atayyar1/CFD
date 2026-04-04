function [Fe, Fw, Fn, Fs] = computeFaceFlux(rho, u, v, Se, Sw, Sn, Ss, Nx, Ny)
% Computes convective face fluxes F_f = rho*v*S_f
Fe= rho .* (u .* Se(:,:,1) + v .* Se(:,:,2));
Fw= rho .* (u .* Sw(:,:,1) + v .* Sw(:,:,2));
Fn= rho .* (u .* Sn(:,:,1) + v .* Sn(:,:,2));
Fs= rho .* (u .* Ss(:,:,1) + v .* Ss(:,:,2));
% I kept it like this, I can simply say that rho=1, but I said lets keep
% like this for any flow of variable density (I know we are not accounting
% for the terms related to the incompressible flows but still.

end