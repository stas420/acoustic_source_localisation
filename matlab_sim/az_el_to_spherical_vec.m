function d_hat = az_el_to_spherical_vec(phi, theta)
    % (azymut, wzniesienie) -> wektory kierunku 3D we wsp. sferycznych
    %
    % Inputs:
    %   phi   - azymut [rad], wektor [N x 1] lub [1 x N]
    %   theta - wzniesienie [rad], wektor [N x 1] lub [1 x N]
    %
    % Output:
    %   d_hat - wektory kierunku [3 x N]
    %           każda KOLUMNA to jeden wektor jednostkowy
    
    phi = phi(:);      % [N x 1]
    theta = theta(:);  % [N x 1]
    N = length(phi);
    
    d_hat = zeros(3, N);
    d_hat(1, :) = cos(theta') .* cos(phi');   % x
    d_hat(2, :) = cos(theta') .* sin(phi');   % y
    d_hat(3, :) = sin(theta');                % z
    
    
    % tak chyba tez jest g:
    % d_hat = [cos(theta) .* cos(phi), ...
    %          cos(theta) .* sin(phi), ...
    %          sin(theta)]';  % [N x 3]' = [3 x N]
end