function [az_norm, el_norm] = normalize_doa_angles(az_est, el_est, az_true, el_true)
    % NORMALIZUJE KĄTY ESTYMOWANE TAK, ABY MINIMALIZOWAĆ BŁĄD WZGLĘDEM PRAWDY
    % Obsługuje przypadki, gdy algorytm zwraca kąty z przeciwnym znakiem
    % ze względu na niejednoznaczność geometryczną
    %
    % Inputs:
    %   az_est, el_est - estymowane kąty [stopnie]
    %   az_true, el_true - prawdziwe kąty [stopnie]
    %
    % Outputs:
    %   az_norm, el_norm - znormalizowane kąty [stopnie]
    %
    % Algorytm:
    %   Sprawdza 4 możliwe kombinacje znaków i wybiera tę,
    %   która minimalizuje odległość kątową na sferze
    
    % Przypadki do sprawdzenia (4 możliwości ze względu na symetrię):
    % 1. (az, el) - bez zmian
    % 2. (-az, el) - odbicie w azymucie
    % 3. (az, -el) - odbicie w elewacji
    % 4. (-az, -el) - odbicie w obu
    
    candidates = [
        az_est,  el_est;
        -az_est, el_est;
        az_est,  -el_est;
        -az_est, -el_est
    ];
    
    % Oblicz błąd dla każdego kandydata
    errors = zeros(4, 1);
    for i = 1:4
        az_cand = candidates(i, 1);
        el_cand = candidates(i, 2);
        
        % Błąd jako odległość kątowa na sferze
        errors(i) = angular_distance(az_true, el_true, az_cand, el_cand);
    end
    
    % Wybierz kandydata z najmniejszym błędem
    [~, best_idx] = min(errors);
    az_norm = candidates(best_idx, 1);
    el_norm = candidates(best_idx, 2);
end

function ang_dist = angular_distance(az1, el1, az2, el2)
    % Oblicza odległość kątową między dwoma punktami na sferze
    % przy użyciu iloczynu skalarnego wektorów jednostkowych
    
    az1_rad = deg2rad(az1);
    el1_rad = deg2rad(el1);
    az2_rad = deg2rad(az2);
    el2_rad = deg2rad(el2);
    
    % Konwersja do wektorów jednostkowych (sferyczne -> kartezjańskie)
    v1 = [cos(el1_rad)*cos(az1_rad); cos(el1_rad)*sin(az1_rad); sin(el1_rad)];
    v2 = [cos(el2_rad)*cos(az2_rad); cos(el2_rad)*sin(az2_rad); sin(el2_rad)];
    
    % Odległość kątowa przez iloczyn skalarny
    % cos(theta) = v1 · v2 / (|v1| * |v2|)
    % Ponieważ v1 i v2 są jednostkowe: cos(theta) = v1 · v2
    cos_angle = dot(v1, v2);
    cos_angle = max(-1, min(1, cos_angle)); % clamp dla stabilności numerycznej
    ang_dist = rad2deg(acos(cos_angle));
end