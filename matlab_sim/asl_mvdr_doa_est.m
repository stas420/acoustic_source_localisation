function [az, el, spec] = asl_mvdr_doa_est(rxSignal, micsPositions, ...
                                        gridAnglePairs, gridSpherical, ...
                                        fs, lambda0)
    %% params setup
    C = 343;
    f_0 = C/lambda0;
    f_bw_half = 250/2;
    [L, M] = size(rxSignal);
    G = size(gridAnglePairs, 1);
    spec = zeros(G, 1);
    epsilon = 1e-6;
       
    %% passband filter to obtain narrowband sig from lambda0 arg
    [b, a] = butter(4, [((f_0 - f_bw_half)/(fs/2)), ...
                        ((f_0 + f_bw_half)/(fs/2))], 'bandpass');
    rxSignal_filt = zeros(size(rxSignal));
    for m = 1:M
        rxSignal_filt(:, m) = filter(b, a, rxSignal(:, m));
    end

    %% autocorr matrix
    Rx = (rxSignal_filt' * rxSignal_filt) / L;
    Rx = Rx + epsilon * eye(M);
    % diagonal loading? co to jest
    Rx_inv = inv(Rx);

    %% MVDR spectrum calc, start with steering vecs
    a = zeros(M, 1);
    for g = 1:G
        for m = 1:M
            tau = micsPositions(:, m)' * gridSpherical(:, g) / C;
            a(m) = exp(1j * 2 * pi * f_0 * tau);
        end

        %a = a / sqrt(M); % normalize?
        spec(g) = 1 / real(a' * Rx_inv * a);
    end

    [~, max_idx] = max(spec(:));
    az = rad2deg(gridAnglePairs(max_idx, 1));
    el = rad2deg(gridAnglePairs(max_idx, 2));

end