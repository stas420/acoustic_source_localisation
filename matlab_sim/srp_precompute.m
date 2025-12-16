function precomp = srp_precompute(micPos, gridSpherical, fMin, fMax, samplingRate, L, C)
    % PRE-OBLICZA STAŁE WARTOŚCI DLA SRP-PHAT
    % Należy wywołać raz przed wieloma uruchomieniami SRP-PHAT
    %
    % Inputs:
    %   micPos - pozycje mikrofonów [3 x M]
    %   gridSpherical - wektory kierunków [3 x G]
    %   fMin, fMax - zakres częstotliwości [Hz]
    %   samplingRate - częstotliwość próbkowania [Hz]
    %   L - długość ramki
    %   C - prędkość dźwięku [m/s]
    %
    % Output:
    %   precomp - struktura z pre-obliczonymi wartościami:
    %       .taus - macierz TDoA [n_pairs x n_dirs]
    %       .all_mic_pairs_idxs - indeksy par mikrofonów [n_pairs x 2]
    %       .freq_indices - indeksy częstotliwości do analizy
    %       .F - wektor częstotliwości
    
    % TDoA matrix - dla każdej pary mikrofonów i każdego kierunku
    all_mic_pairs_idxs = nchoosek(1:size(micPos, 2), 2);
    n_pairs = size(all_mic_pairs_idxs, 1);
    n_dirs = size(gridSpherical, 2);
    taus = zeros(n_pairs, n_dirs);
    
    for p = 1:n_pairs
        l = all_mic_pairs_idxs(p, 1);
        m = all_mic_pairs_idxs(p, 2);
        delta_v = micPos(:, m) - micPos(:, l);  % [3 x 1]
        % delta_v' * gridSpherical = [1 x 3] * [3 x n_dirs] = [1 x n_dirs]
        taus(p, :) = (delta_v' * gridSpherical) / C;
    end
    
    % Freq bins setup - które biny FFT analizować
    freqs = (0:L-1) * (samplingRate / L);  % [L x 1]
    freq_mask = (freqs >= fMin) & (freqs <= fMax);
    freq_indices = find(freq_mask);
    F = freqs(freq_indices);
    
    % Zwróć strukturę z pre-obliczonymi wartościami
    precomp.taus = taus;
    precomp.all_mic_pairs_idxs = all_mic_pairs_idxs;
    precomp.freq_indices = freq_indices;
    precomp.F = F;
    
    fprintf('SRP-PHAT pre-compute: %d par mikrofonów, %d kierunków, %d częstotliwości\n', ...
        n_pairs, n_dirs, length(F));
end