function [M, M_error, p, p_error] = supersonic_theo(P, P_error, y)
     %{
     Input
     pressure : vector(2,1)
        Static and total pressures corresponding to port 1.
        (1,:): Static port pressures
          (2,:): Total pressures
     pressure_err : vector (2, 7)
        Precision uncertainty of static and total pressures.
        (1,:) : Static port pressures errors
        (2,:) : Total pressures errors
     y : float
        Ratio of specific heats (gamma).
    
    Returns
    M : matrix (1, 7)
        The Mach numbers for the input pressure data.
    M_error : matrix (1, 7)
        The Mach number error
    p : matrix (1, 7)
        The pressure for the theoretical data.
    p_error : matrix (1, 7)
        error associated with all pressure values at each port
    %}
    Ar = [1.06, 1.00, 1.05, 1.15, 1.23, 1.27, 1.28];
    
    p_error = zeros(1, 7);
    p = zeros(1, 7);
    M = zeros(1,8);
    M_error = zeros(1,8);

    for i = 1:size(Ar)
        f = @(m) 1/m *(((2/(y+1)*(1 + ((y-1)/2)*m^2))^((y+1)/2*(y-1)))) - Ar(i);
        if i == 1
            % since flow will likely be subsonic before throat
            [M(i), M_error(i)] = fzero(f, 0.8);
        elseif i == 2
            % since flow will likely be just over Mach 1 at the second tap
            % since its the throat
            [M(i), M_error(i)] = fzero(f, 1.1);
        else 
            % Supersonic everywhere else
            [M(i), M_error(i)] = fzero(f, 1.2);
        end
    end



    g = @(m) (((y+1).*m.^2./2).^(-y./(y-1))).*...
             (((2.*y.*m.^2)./(y+1))-(y-1)./(y+1)).^(1./y-1);
    p(1, :) = P(2,1)*g(M(1,:));


    dp = @(p,m) p.*(-y.*m.^((-3.*y+1)./y-1).*(2.*y.*m.^2 - y + 1).*(y+1).^...
            ((-y-1)./y-1)) + (2.*y.*m.^((-y-1)/y-1) .* (2.*y.*m.^2 - y+1)...
            .^((-y+2)./y-1) .*(y+1).^((-y-1)./y-1)./(2.^((-y+1)./(y-1).*(y-1))));

    p_error = sqrt((P_error(2).*g(M)).^2 + (M_error.*dp(P, M)).^2);
end 

        

