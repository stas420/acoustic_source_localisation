function [az, el, info] = asl_das_doa_est(rxSignal, micsPositions, ...
                            gridAnglePairs, gridDirVectors, fs)
    % POPRAWIONA WERSJA - naprawiono obliczenia tau
    
    C = 343;
    freq_min = 300;
    freq_max = 3300;

    G = size(gridAnglePairs, 1);
    [L, M] = size(rxSignal);
    fft_bin_min = floor(freq_min * L / fs);
    fft_bin_max = floor(freq_max * L / fs);
    fft_bins_n = fft_bin_max - fft_bin_min + 1;
    rx_fft = zeros(fft_bins_n, M);
    w = 1/M;
    taus = zeros(G, M);

    % POPRAWKA: pierwszy mikrofon jako referencja
    % dla każdego kierunku i każdego mikrofonu obliczamy opóźnienie
    for m = 1:M
        % Różnica pozycji względem pierwszego mikrofonu
        d = (micsPositions(:, m) - micsPositions(:, 1));  % [3 x 1]
        
        % Dla każdego kierunku w siatce
        % gridDirVectors to [3 x G], więc gridDirVectors' to [G x 3]
        % d to [3 x 1]
        % Wynik: [G x 3] * [3 x 1] = [G x 1]
        taus(:, m) = (gridDirVectors' * d) / C;  % FIX: poprawne mnożenie macierzowe
    end
    
    % Pierwszy mikrofon ma oczywiście tau = 0 dla wszystkich kierunków
    % (to już jest zapewnione powyżej, bo d = 0 dla m=1)

    % FFT każdego kanału
    for m = 1:M
        tmp_fft = fft(rxSignal(:, m));
        rx_fft(:,m) = tmp_fft(fft_bin_min:fft_bin_max);
    end

    % Dla każdego kierunku w siatce
    for g = 1:G
        freq_spec = zeros(fft_bins_n, 1);
        freqs = (fft_bin_min:fft_bin_max)' * fs / L;  % częstotliwości dla tych binów
        
        for m = 1:M
            % Delay-and-sum: kompensujemy opóźnienie dla danego kierunku
            freq_spec = freq_spec + w * rx_fft(:, m) .* ...
                        exp(-1j * 2 * pi * freqs * taus(g, m));
        end
        
        % Moc sygnału dla tego kierunku
        power_spec(g) = sum(abs(freq_spec).^2);
    end

    % Znajdź maksimum
    [~, max_idx] = max(power_spec);
    az = rad2deg(gridAnglePairs(max_idx, 1));
    el = rad2deg(gridAnglePairs(max_idx, 2));
    
    % Zwróć dodatkowe info
    info = struct('maxPower', power_spec(max_idx), ...
                  'maxIndex', max_idx, ...
                  'powerSpectrum', power_spec);  % dodane pełne spektrum
end