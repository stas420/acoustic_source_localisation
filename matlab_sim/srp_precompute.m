function precomp = srp_precompute(micPos, gridSpherical, fMin, fMax, samplingRate, L, C)   
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
    
    % Freq bins setup
    freqs = (0:L-1) * (samplingRate / L);  % [L x 1]
    freq_mask = (freqs >= fMin) & (freqs <= fMax);
    freq_indices = find(freq_mask);
    F = freqs(freq_indices);
    
    % ret
    precomp.taus = taus;
    precomp.all_mic_pairs_idxs = all_mic_pairs_idxs;
    precomp.freq_indices = freq_indices;
    precomp.F = F;

end