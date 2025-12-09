function [az, el, info] = asl_das_doa_est(rxSignal, micsPositions, ...
                            gridAnglePairs, gridDirVectors, fs)

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
    power_spec = zeros(G, 1);

    % first mic here is left 0 as a reference one
    for m = 2:M
        d = (micsPositions(:, m) - micsPositions(:, 1));
        taus(:, m) = gridDirVectors(:, :)' * ... 
                    d / C; 
    end

    for m = 1:M
        tmp_fft = fft(rxSignal(:, m));
        rx_fft(:,m) = tmp_fft(fft_bin_min:fft_bin_max);
    end

    for g = 1:G
        freq_spec = zeros(fft_bins_n, 1);
        for m = 1:M
            freqs = (fft_bin_min:fft_bin_max)' * fs / L;
            freq_spec(:) = freq_spec(:) +  w * rx_fft(:, m) ...
                            .* exp(-1j * 2 * pi * freqs * taus(g,m));
        end
        power_spec(g) = sum(abs(freq_spec).^2);
    end

    [~, max_idx] = max(power_spec);
    az = rad2deg(gridAnglePairs(max_idx, 1));
    el = rad2deg(gridAnglePairs(max_idx, 2));
    info = struct('maxPower', power_spec(max_idx), 'maxIndex', max_idx);
end