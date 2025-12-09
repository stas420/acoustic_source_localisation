function [az, el, info] = asl_srp_phat(rxSignal, micPos, fMin, fMax, ...
                                       samplingRate, gridAnglePairs, ...
                                       gridSpherical, C, L)
    %% --- TDoA matrix --- << this ood be considered precalc in other situations
    % ---- for every d vector in D calculate TDoA per every possible mic pair
    all_mic_pairs_idxs = nchoosek(1:size(micPos, 2), 2); % Generate all mic pair indices
    n_pairs = size(all_mic_pairs_idxs, 1); % Number of mic pairs
    n_dirs = size(gridSpherical, 2); % Number of direction vectors
    taus = zeros(n_pairs, n_dirs);
    
    for p = 1:n_pairs
        l = all_mic_pairs_idxs(p, 1);
        m = all_mic_pairs_idxs(p, 2);
        
        % difference of positions in cartesian coords
        delta_v = micPos(:, m) - micPos(:, l);  % --> [3 x 1]
        
        % delta_v' * D = [1 x 3] * [3 x n_dirs] = [1 x n_dirs]
        taus(p, :) = (delta_v' * gridSpherical) / C; % <-- this is from that weird formula considering direction angles only
    end
    
    %fprintf('TDoA matrix: [%d x %d]\n', n_pairs, n_dirs);

    %% --- freq bins setup ---  << this too
    % FFT's bins frequencies 
    freqs = (0:L-1) * (samplingRate / L);  % [L x 1]
    % and here we keep only the band of interest
    freq_mask = (freqs >= fMin) & (freqs <= fMax);
    freq_indices = find(freq_mask);

    % fprintf("fft after setup, starting\n");

    %% mic frames' DFT and initial filtering << and this also!
    micFramesFFT = fft(rxSignal, [], 1);  % [L x M]    
    micFramesFFT_filtered = micFramesFFT(freq_indices, :);  % [n_freqs x M]
    F = freqs(freq_indices);
    n_freqs = length(F);
    
    %fprintf('FFT done: %d frequencies (%.0f - %.0f Hz)\n', ...
    %        n_freqs, F(1), F(end));

    %% --- GCC-PHAT computation ---    
    % GGC-PHAT is a function of frequency per every mic pair considered
    gcc_phat_all = zeros(n_pairs, n_freqs);  % [n_pairs x n_freqs]
    
    %fprintf('Computing GCC-PHAT for all pairs...\n');
    tic
    for p = 1:n_pairs
        l = all_mic_pairs_idxs(p, 1);
        m = all_mic_pairs_idxs(p, 2);
        
        X_l = micFramesFFT_filtered(:, l);  % [n_freqs x 1]
        X_m = micFramesFFT_filtered(:, m);  % [n_freqs x 1]
        
        gcc_phat_all(p, :) = asl_gcc_phat(X_l, X_m);  % [1 x n_freqs]
    end

    %fprintf('GCC-PHAT computed for all pairs\n');

    %% --- SRP map computation ---
    % map is sum of GCCs across every mic pair and freq per every direction
    srp_map = zeros(n_dirs, 1);
    
    %fprintf('Computing SRP map...\n');

    % for every direction considered...
    for u = 1:n_dirs            
        srp_value = 0;
        % ...and for every mic pair...
        for p = 1:n_pairs
            % ...for TDoA value determined by the specific mic pair and
            % direction...
            tau = taus(p, u);
            
            % ...and for each frequency bin considered...
            for f_idx = 1:n_freqs
                f = F(f_idx);
                
                % phase shift for this specific frequency and TDoA
                phase_shift = exp(1j * 2 * pi * f * tau);
                
                % GCC-PHAT value for this particular pair and freq
                g = gcc_phat_all(p, f_idx);
                
                % SRP-PHAT(u, f, xl, xm) = 
                % Re[GCC-PHAT(f,xl,xm) * e^(j2pif * tau_lm(u))
                contrib = real(g * phase_shift);
                % and it's a sum!
                srp_value = srp_value + contrib;
            end
        end
        
        % normalisation per mics pair number - some say do, some say don't...
        srp_map(u) = srp_value; %/n_pairs;
    end
    
    %fprintf('SRP map computed\n');

    %% --- find maximum in SRP map ---
    [~, max_idx] = max(srp_map);
    info.time = toc;
    %fprintf("SRP-PHAT is done\n");
    az = rad2deg(gridAnglePairs(max_idx, 1));
    el = rad2deg(gridAnglePairs(max_idx, 2));
end