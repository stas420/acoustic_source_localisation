function [az, el, info] = asl_srp_phat(rxSignal, micPos, fMin, fMax, ...
                                       samplingRate, gridAnglePairs, ...
                                       gridSpherical, C, L, precomp)
    % SRP-PHAT z opcjonalnymi pre-computacjami dla przyspieszenia
    % 
    % Inputs:
    %   rxSignal - sygnały z mikrofonów [L x M]
    %   micPos - pozycje mikrofonów [3 x M]
    %   fMin, fMax - zakres częstotliwości [Hz]
    %   samplingRate - częstotliwość próbkowania [Hz]
    %   gridAnglePairs - siatka kątów [az, el] w radianach [G x 2]
    %   gridSpherical - wektory kierunków [3 x G]
    %   C - prędkość dźwięku [m/s]
    %   L - długość ramki
    %   precomp (opcjonalny) - struktura z pre-obliczonymi wartościami
    %
    % Outputs:
    %   az, el - estymowane kąty [stopnie]
    %   info - struktura z dodatkowymi informacjami
    
    %% --- Pre-computacje lub ich użycie ---
    if nargin < 10 || isempty(precomp)
        % Oblicz TDoA matrix (jeśli nie podano pre-compute)
        all_mic_pairs_idxs = nchoosek(1:size(micPos, 2), 2);
        n_pairs = size(all_mic_pairs_idxs, 1);
        n_dirs = size(gridSpherical, 2);
        taus = zeros(n_pairs, n_dirs);
        
        for p = 1:n_pairs
            l = all_mic_pairs_idxs(p, 1);
            m = all_mic_pairs_idxs(p, 2);
            delta_v = micPos(:, m) - micPos(:, l);
            taus(p, :) = (delta_v' * gridSpherical) / C;
        end
        
        % Freq bins setup
        freqs = (0:L-1) * (samplingRate / L);
        freq_mask = (freqs >= fMin) & (freqs <= fMax);
        freq_indices = find(freq_mask);
        F = freqs(freq_indices);
    else
        % Użyj pre-obliczonych wartości (szybsze!)
        taus = precomp.taus;
        all_mic_pairs_idxs = precomp.all_mic_pairs_idxs;
        freq_indices = precomp.freq_indices;
        F = precomp.F;
        n_pairs = size(all_mic_pairs_idxs, 1);
        n_dirs = size(gridSpherical, 2);
    end
    
    n_freqs = length(F);

    %% --- DFT i filtrowanie ---
    micFramesFFT = fft(rxSignal, [], 1);  % [L x M]
    micFramesFFT_filtered = micFramesFFT(freq_indices, :);  % [n_freqs x M]

    %% --- GCC-PHAT dla wszystkich par ---
    gcc_phat_all = zeros(n_pairs, n_freqs);  % [n_pairs x n_freqs]
    
    for p = 1:n_pairs
        l = all_mic_pairs_idxs(p, 1);
        m = all_mic_pairs_idxs(p, 2);
        
        X_l = micFramesFFT_filtered(:, l);  % [n_freqs x 1]
        X_m = micFramesFFT_filtered(:, m);  % [n_freqs x 1]
        
        gcc_phat_all(p, :) = asl_gcc_phat(X_l, X_m);  % [1 x n_freqs]
    end

    %% --- SRP map - WEKTORYZACJA (zamiast 3 pętli tylko 1!) ---
    % Stara wersja miała:
    %   for u = 1:n_dirs
    %       for p = 1:n_pairs
    %           for f_idx = 1:n_freqs
    % Nowa wersja:
    %   for u = 1:n_dirs
    %       [operacje macierzowe]
    
    srp_map = zeros(n_dirs, 1);
    
    % Dla każdego kierunku oblicz SRP wartość
    for u = 1:n_dirs
        % Phase shifts dla wszystkich par i częstotliwości naraz
        % taus(p, u) - opóźnienia dla kierunku u i pary p
        % F to [n_freqs x 1] - częstotliwości
        % Wynik: [n_pairs x 1] * [1 x n_freqs] = [n_pairs x n_freqs]
        
        tau_u = taus(:, u);  % [n_pairs x 1]
        tau_u = tau_u(:);    % zapewnij że to kolumna [n_pairs x 1]
        F_row = F(:)';       % zapewnij że to wiersz [1 x n_freqs]
        % Broadcasting: [n_pairs x 1] * [1 x n_freqs] = [n_pairs x n_freqs]
        phase_shifts = exp(1j * 2 * pi * tau_u * F_row);  % [n_pairs x n_freqs]
        
        % Element-wise multiplication i suma
        % gcc_phat_all: [n_pairs x n_freqs]
        % phase_shifts: [n_pairs x n_freqs]
        contrib = real(gcc_phat_all .* phase_shifts);  % [n_pairs x n_freqs]
        
        % Suma po wszystkich parach i częstotliwościach
        srp_map(u) = sum(contrib(:));
    end

    %% --- Znajdź maksimum ---
    [~, max_idx] = max(srp_map);
    az = rad2deg(gridAnglePairs(max_idx, 1));
    el = rad2deg(gridAnglePairs(max_idx, 2));
    
    % Dodatkowe info
    info.srp_map = srp_map;
    info.max_idx = max_idx;
end