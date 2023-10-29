function [M, M_error, p, p_error] = subsonic_theo(P, P_error, y)
    %{
    Inputs
    P : matrix (2, 7)
        Static and total pressures corresponding to each position.
        (1,:): Static port pressures
        (2,:): Total pressures
    P_error : matrix (2, 7)
        Precision uncertainty of static and total pressures
        (1,:): Static port pressures
        (2,:): Total pressures
    y : float
        Ratio of specific heats (gamma).

    Returns
    M : matrix (1, 7)
        The Mach numbers for the input pressure data.
    M_error : matrix (1, 7)
        The Mach number uncertainties.
    p : matrix (1, 7)
        The pressure for the theoretical data.
    p_err : matrix(1, 7)
        Uncertainties associated with the pressures.
    %}

    Ar = [1.06, 1.00, 1.05, 1.15, 1.23, 1.27, 1.28];
    
    p_error = zeros(1, 7);
    p = zeros(1, 7);
    M = zeros(1,7);
    M_error = zeros(1,7);

    for i = 1:size(Ar)
        f = @(m) 1/m *(((2/(y+1)*(1 + ((y-1)/2)*m^2))^((y+1)/2*(y-1)))) - Ar(i);
        if i == 1
            % since flow will likely be subsonic before throat
            [M(i), M_error(i)] = fzero(f, 0.8);
        elseif i == 2
            % since flow will likely be just over Mach 1 at the second tap
            % since its the throat
            [M(i), M_error(i)] = fzero(f, 1);
        else 
            % Subsonic everywhere else
            [M(i), M_error(i)] = fzero(f, 0.4);
        end
    end
    
    g = @(m)(1 + (y - 1)./2 .*m.^2).^(-y ./(y - 1));
    
    p(1, :) = P(2,1) * g(M(1,:));
    
    dp = @(m)-(2.*m.*y*(y./2 - 1./2))./(((y./2 - 1./2).*m.^2 + 1) ...
        .^(y./(y - 1) + 1).*(y - 1));

    p_error = sqrt((g(M) .* P_error(2)).^2 + (dp(M) .* M_error).^2);

    
end
