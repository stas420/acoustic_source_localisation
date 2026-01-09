function [az, el, info] = asl_das_doa_est(rxSignal, gridAnglePairs, ...
                                          fs, taus)
    %% -- info
    % Delay-and-Sum beamforming used for DoA estimation is simply
    % just grid search by applying proper time delays per direction
    % and mic. In time domain, it is quite tough, due to interpolation
    % operations needed. Generally it is done through FFT and magnitude
    % computation (optionally IFFT, but rarely) because it's simpler 
    % both computationally and conceptually.
    
    %% -- setup
    C = 343;
    G = size(gridAnglePairs, 1);
    [L, M] = size(rxSignal);
    %w = 1/(M*M);
    
    %% -- compute fft
    rx_fft = fft(rxSignal);

    %% -- compute power per direction with bandpass filtering
    f_min = 300;
    f_max = 3000;
    freqs_all = (0:(L-1))' * fs/L;
    freq_mask = (freqs_all >= f_min) & (freqs_all <= f_max);
    freqs_filt = freqs_all(freq_mask);
    rx_fft_filt = rx_fft(freq_mask, :);
    
    power_spec = zeros(G, 1);

    % for each direction g...
    for g = 1:G
        Y_sum = zeros(length(freqs_filt), 1);

        for m = 1:M
            Y_sum = Y_sum + ...
                rx_fft_filt(:, m) .* exp(-1j * 2 * pi * freqs_filt * taus(g, m));
        end
        
        power_spec(g) = sum(real(Y_sum));

    end

    power_spec = power_spec / max(abs(power_spec));

    %% -- find max 
    [~, max_idx] = max(power_spec);
    az = rad2deg(gridAnglePairs(max_idx, 1));
    el = rad2deg(gridAnglePairs(max_idx, 2));
    
    %% -- other stuff
    info = struct('maxPower', power_spec(max_idx), ...
                  'maxIndex', max_idx, ...
                  'powerSpectrum', power_spec);
end