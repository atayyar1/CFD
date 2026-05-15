function [ge,gw,gn,gs] = getInterpolationFactors(Nx,Ny,Vc,CE,CW,CN,CS,method)

% Initialize arrays
ge = zeros(Nx,Ny);
gw = zeros(Nx,Ny);
gn = zeros(Nx,Ny);
gs = zeros(Nx,Ny);

if strcmp(method,'volume')

    % Volume-based interpolation
    for i = 1:Nx
    for j = 1:Ny

        if i < Nx
            ge(i,j) = Vc(i,j)/(Vc(i,j) + Vc(i+1,j));
        end
        if i > 1
            gw(i,j) = Vc(i,j)/(Vc(i,j) + Vc(i-1,j));
        end
        if j < Ny
            gn(i,j) = Vc(i,j)/(Vc(i,j) + Vc(i,j+1));
        end

        if j > 1
            gs(i,j) = Vc(i,j)/(Vc(i,j) + Vc(i,j-1));
        end

    end
    end

elseif strcmp(method,'distance')

    % Distance-based interpolation
    for i = 1:Nx
    for j = 1:Ny

        if i < Nx
            dP = norm(squeeze(CE(i,j,:)));
            dE = norm(squeeze(CW(i+1,j,:)));
            ge(i,j) = dE/(dP + dE);
        end

        if i > 1
            dP = norm(squeeze(CW(i,j,:)));
            dW = norm(squeeze(CE(i-1,j,:)));
            gw(i,j) = dW/(dP + dW);
        end

        if j < Ny
            dP = norm(squeeze(CN(i,j,:)));
            dN = norm(squeeze(CS(i,j+1,:)));
            gn(i,j) = dN/(dP + dN);
        end

        if j > 1
            dP = norm(squeeze(CS(i,j,:)));
            dS = norm(squeeze(CN(i,j-1,:)));
            gs(i,j) = dS/(dP + dS);
        end

    end
    end

else
    error('Please choose the interpolation scheme, "distance" or "volume"')

end

end
%% CHat gpt suggested that I use squeeze to fix the fact that Matlab is storing the vectors as higher dimension
%% Also, the strmp choice was suggested by chat gpt