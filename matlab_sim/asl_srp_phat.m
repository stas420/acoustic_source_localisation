function [az, el, info] = asl_srp_phat(rxSignal, micPos, fMin, fMax, ...
                                       samplingRate, gridAnglePairs, ...
                                       gridSpherical, C, L, precomp)
    %% --- LUT ---
    if nargin < 10 || isempty(precomp)
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
        
        freqs = (0:L-1) * (samplingRate / L);
        freq_mask = (freqs >= fMin) & (freqs <= fMax);
        freq_indices = find(freq_mask);
        F = freqs(freq_indices);
    else
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

    %% --- SRP map  ---    
    srp_map = zeros(n_dirs, 1);
    
    for u = 1:n_dirs
        
        tau_u = taus(:, u);  % [n_pairs x 1]
        tau_u = tau_u(:);
        F_row = F(:)';

        phase_shifts = exp(1j * 2 * pi * tau_u * F_row);  % [n_pairs x n_freqs]
        contrib = real(gcc_phat_all .* phase_shifts);  % [n_pairs x n_freqs]
        srp_map(u) = sum(contrib(:)) / (n_pairs * n_freqs);
    end

    %% --- maksimum ---
    [~, max_idx] = max(srp_map);
    az = rad2deg(gridAnglePairs(max_idx, 1));
    el = rad2deg(gridAnglePairs(max_idx, 2));
    
    info.srp_map = srp_map;
    info.max_idx = max_idx;
end