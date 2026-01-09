function [az, el, info] = asl_music_2d(rxSignal, lambda0, micPositions, ...
                                       numSources, gridOfAnglesRad, ...
                                       gridSpherical, fs)
    M = size(micPositions, 2);
    
    if size(rxSignal, 1) ~= M
        rxSignal = rxSignal';
    end
    
    %% Filtrowanie
    C = 343;
    fc = C/lambda0;
    bandwidth_hz = 500;
    f_low = fc - bandwidth_hz/2;
    f_high = fc + bandwidth_hz/2;
    
    rxSig_fft = fft(rxSignal, [], 2);
    freq_axis = (0:size(rxSignal, 2)-1) * fs / size(rxSignal, 2);
    freq_mask = (freq_axis >= f_low) & (freq_axis <= f_high);
    rxSig_fft(:, ~freq_mask) = 0;
    rxSignal_filtered = ifft(rxSig_fft, [], 2);  % USUŃ real()!
    
    %% Eigen decomp
    Rxx = (rxSignal_filtered * rxSignal_filtered') / size(rxSignal_filtered, 2);
    Rxx = (Rxx + Rxx') / 2;
    [V, D] = eig(Rxx);
    [~, idx] = sort(diag(D), 'descend');
    V = V(:, idx);
    U_N = V(:, (numSources+1):end);
    
    eigenvals = diag(D);
    eigenvals = eigenvals(idx);
    
    %% MUSIC spectrum
    f_0 = C / lambda0;  % częstotliwość centralna [Hz]
    n_dirs = size(gridOfAnglesRad, 1);
    P_MUSIC = zeros(n_dirs, 1);
    U_N_proj = U_N * U_N';
    
    mic_centre = mean(micPositions, 2);
    
    for n = 1:n_dirs
        dir_vec = gridSpherical(:, n);
        
        % Steering vector - POPRAWIONE!
        a = zeros(M, 1);
        for m = 1:M
            delta_pos = micPositions(:, m) - mic_centre;
            tau = (delta_pos' * dir_vec) / C;
            a(m) = exp(1j * 2 * pi * f_0 * tau);  % 2πf₀τ zamiast k*τ
        end
        
        a = a / sqrt(a' * a);
        P_MUSIC(n) = 1 / real(a' * U_N_proj * a);
    end
    
    %% Finding maxes
    [~, max_idx] = maxk(P_MUSIC, min(numSources+2, length(P_MUSIC)));
    angles = gridOfAnglesRad(max_idx(1:numSources), :);
    az = rad2deg(angles(:, 1));
    el = rad2deg(angles(:, 2));
    
    info.P_MUSIC = P_MUSIC;
    info.eigenvalues = eigenvals;
end