function [az, el, spec] = asl_mvdr_doa_est(rxSignal, micsPositions, ...
                                        gridAnglePairs, gridSpherical, ...
                                        fs, lambda0)
    C = 343;
    f_0 = C/lambda0;
    f_bw_half = 125;
    f_low = f_0 - f_bw_half;
    f_high = f_0 + f_bw_half;

    [N_samples, M] = size(rxSignal);  % N_samples=1024, M=16 mikrofonów
    G = size(gridAnglePairs, 1);
    spec = zeros(G, 1);
    
    %% Filtrowanie - ZOSTAW ZESPOLONE!
    rxSig_fft = fft(rxSignal, [], 1);  % FFT wzdłuż sampli (dim=1)
    freq_axis = (0:N_samples-1) * fs / N_samples;
    freq_mask = (freq_axis >= f_low) & (freq_axis <= f_high);
    rxSig_fft(~freq_mask, :) = 0;  % Zeruj wiersze (częstotliwości)
    rxSignal_filt = ifft(rxSig_fft, [], 1);  % NIE bierz real()!

    %% Macierz korelacji [M × M]
    Rx = (rxSignal_filt' * rxSignal_filt) / N_samples;
    % rxSignal_filt: [1024 × 16]
    % rxSignal_filt': [16 × 1024] (transpozycja sprzężona)
    % Iloczyn: [16 × 1024] × [1024 × 16] = [16 × 16] ✓
    
    Rx = (Rx + Rx') / 2;  % Hermitianizacja
    
    epsilon = 1e-4 * trace(real(Rx)) / M;
    Rx = Rx + epsilon * eye(M);
    Rx_inv = inv(Rx);

    %% Steering vectors
    mic_centre = mean(micsPositions, 2);
    
    for g = 1:G
        dir_vec = gridSpherical(:, g);
        
        a = zeros(M, 1);
        for m = 1:M
            delta_pos = micsPositions(:, m) - mic_centre;
            tau = (delta_pos' * dir_vec) / C;
            a(m) = exp(-1j * 2 * pi * f_0 * tau);
        end
        
        a = a / sqrt(a' * a);
        spec(g) = 1 / real(a' * Rx_inv * a);
    end

    % [~, max_idx] = max(spec(:));
    [~, XI] = maxk(spec(:), 2);
    max_idx = XI(1);
    az = rad2deg(gridAnglePairs(max_idx, 1));
    el = rad2deg(gridAnglePairs(max_idx, 2));
end