function bDC = computeSMARTCorrection(phi, Fe, Fw, Fn, Fs,dphidx, dphidy, XC, YC, Xe, Ye, Xw, Yw, Xn, Yn, Xs, Ys, Nx, Ny)
bDC = zeros(Nx, Ny);

for i = 1:Nx
    for j = 1:Ny
        % --- east face ---
        if i < Nx
 if Fe(i,j) > 0
 % flow goes left to right, C is (i,j), D is (i+1,j)
                phiC = phi(i,j);
                phiD = phi(i+1,j);
                % phiU from professor's formula: phiD - grad_C * 2*r_C
                      dx = XC(i+1,j) - XC(i,j);
                dy = YC(i+1,j) - YC(i,j);
                phiU = phiD - 2*(dphidx(i,j)*dx + dphidy(i,j)*dy);
            else
                % flow goes right to left, C is (i+1,j), D is (i,j)
                phiC = phi(i+1,j);
                phiD = phi(i,j);
                dx = XC(i,j) - XC(i+1,j);
                dy = YC(i,j) - YC(i+1,j);
                phiU = phiD - 2*(dphidx(i+1,j)*dx + dphidy(i+1,j)*dy);
            end

            phi_smart = SMARTlimiter(phiC, phiD, phiU);
            phi_uw = phiC; % upwind is just the C cell value

            % conservative: what leaves one cell enters the other
            bDC(i,j)   = bDC(i,j)   - Fe(i,j)*(phi_smart - phi_uw);
            bDC(i+1,j) = bDC(i+1,j) + Fe(i,j)*(phi_smart - phi_uw);
        end

        % --- north face ---
        if j < Ny
            if Fn(i,j) > 0
                phiC = phi(i,j);
                phiD = phi(i,j+1);
                dx = XC(i,j+1) - XC(i,j);
                dy = YC(i,j+1) - YC(i,j);
                phiU = phiD - 2*(dphidx(i,j)*dx + dphidy(i,j)*dy);
            else
                phiC = phi(i,j+1);
                phiD = phi(i,j);
                dx = XC(i,j) - XC(i,j+1);
                dy = YC(i,j) - YC(i,j+1);
                phiU = phiD - 2*(dphidx(i,j+1)*dx + dphidy(i,j+1)*dy);
            end

            phi_smart = SMARTlimiter(phiC, phiD, phiU);
            phi_uw = phiC;

            bDC(i,j)   = bDC(i,j)   - Fn(i,j)*(phi_smart - phi_uw);
            bDC(i,j+1) = bDC(i,j+1) + Fn(i,j)*(phi_smart - phi_uw);
        end

    end
end

end

function phi_f = SMARTlimiter(phiC, phiD, phiU)

denom = phiD - phiU;
if abs(denom) < 1e-10
    phi_f = phiC;
    return
end

phiC_t = (phiC - phiU) / denom;

if phiC_t < 0 || phiC_t > 1
    phi_f = phiC;                               % upwind

elseif phiC_t < 1/6
    phi_f = phiU + denom * 3*phiC_t;           % slope 3 region

elseif phiC_t <= 7/10
    phi_f = phiU + denom*(3/4*phiC_t + 3/8);  % QUICK region

else
    phi_f = phiU + denom*(1/3*phiC_t + 2/3);  % modified upper region
end

end