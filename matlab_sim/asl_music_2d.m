function [az, el, info] = asl_music_2d(rxSignal, lambda0, micSpacingYZ, ...
                                       numMicsInRowCol, numSources, gridOfAnglesRad)
   
    M = numMicsInRowCol * numMicsInRowCol;

    if size(rxSignal, 1) ~= M
        rxSignal = rxSignal';
    end
    
    % [M, L] = size(rxSignal);
    % bandpass - classic impl of MUSIC is narrowband!
    C = 343;
    fc = C/lambda0;
    fs = 16000;
    bandwidth_hz = fs/100;
    
    f_low = fc - bandwidth_hz/2;
    f_high = fc + bandwidth_hz/2;
    
    if f_low <= 0 || f_high >= fs/2
        fprintf('freqs oob');
        fprintf('fc=%.1f Hz, bw=%.1f Hz, fs=%.1f Hz, nyq=%.1f Hz\n', ...
                fc, bandwidth_hz, fs, fs/2);
        rxSignal_filtered = rxSignal; % << idk what to do, this shouldn't happen
    else
        % Filtracja każdego kanału
        M_elements = size(rxSignal, 1);
        rxSignal_filtered = zeros(size(rxSignal));

        for m = 1:M_elements
            % bandpass(signal, [f_low f_high], fs)
            rxSignal_filtered(m, :) = bandpass(rxSignal(m, :), [f_low f_high], fs);
        end
    end

    rxSignal_ol = rxSignal;
    rxSignal = rxSignal_filtered;

    % f_axis = fs/L*(0:(L-1));
    % figure;
    % subplot(1,2,1);
    % plot(f_axis, fft(rxSignal(1,:)));
    % title('rx after');
    % subplot(1,2,2);
    % plot(f_axis, fft(rxSignal_ol(1,:)));
    % title('rx bef');

    tic
    % autocorr of incoming sig Rxx
    Rxx = (rxSignal * rxSignal')/size(rxSignal, 2);
    Rxx = (Rxx + Rxx')/2;
    
    % eigen decomp
    [V, D] = eig(Rxx); % Eigen decomposition of the autocorrelation matrix
    [~, idx] = sort(diag(D), 'descend'); % Sort eigenvalues in descending order
    V = V(:, idx); % Sort eigenvectors according to sorted eigenvalues
    U_N = V(:, (numSources+1):end); % noise subspace

    eigenvals = diag(D);
    eigenvals = eigenvals(idx);
    % fprintf('eigvals: ');
    % fprintf('%.2e ', eigenvals(1:min(5,length(eigenvals))));
    % fprintf('\n');

    % calc MUSIC pseudo-spectrum
    k = 2 * pi / lambda0;
    n_dirs = size(gridOfAnglesRad, 1);
    P_MUSIC = zeros(n_dirs, 1);
    U_N_proj = U_N * U_N';

    for n = 1:n_dirs
        az_n = gridOfAnglesRad(n, 1);
        el_n = gridOfAnglesRad(n, 2);

        % URA in yz-plane, faced with +x
        % az { +X : +Y }, el { XY-plane up to +Z}
        
        % along y col : sin(az) * cos(el)
        a_y = exp(1j * k * micSpacingYZ * (0:(numMicsInRowCol-1))' * sin(az_n) * cos(el_n));
        
        % along z rows: sin(el)
        a_z = exp(1j * k * micSpacingYZ * (0:(numMicsInRowCol-1))' * sin(el_n));
        
        % kronecker - this depends on mics order << to be reviewed
        a_az_el_n = kron(a_z, a_y);
        
        P_MUSIC(n) = 1 / real(a_az_el_n' * U_N_proj * a_az_el_n);
    end
    
    % finding maxes
    [peak_vals, max_idx] = maxk(P_MUSIC, min(numSources+2, length(P_MUSIC)));
    
    % fprintf('\ntop %d peaks:\n', length(max_idx));
    % for i = 1:length(max_idx)
    %     fprintf('  %d: az = %.1f, el = %.1f, val = %.2e\n', ...
    %             i, rad2deg(gridOfAnglesRad(max_idx(i),1)), ...
    %             rad2deg(gridOfAnglesRad(max_idx(i),2)), peak_vals(i));
    % end
    
    angles = gridOfAnglesRad(max_idx(1:numSources), :);
    az = rad2deg(angles(:, 1));
    el = rad2deg(angles(:, 2));
    
    for i = 1:numSources
        fprintf('  doa %d: az (deg) = %.1f, el (deg) = %.1f\n', i, az(i), el(i));
    end
    
    info.P_MUSIC = P_MUSIC;
    info.eigenvalues = eigenvals;
    
    % visualize_music_spectrum(gridOfAnglesRad, P_MUSIC, az, el, 45, 26.3);
end

% function visualize_music_spectrum(grid, P_MUSIC, detected_az, detected_el, true_az, true_el)
%     az_unique = uniquetol(grid(:,1), 1e-9);
%     el_unique = uniquetol(grid(:,2), 1e-9);
% 
%     if length(az_unique) * length(el_unique) == size(grid, 1)
%         P_2D = reshape(P_MUSIC, length(el_unique), length(az_unique));
%         P_2D_dB = 10*log10(P_2D);
%         P_2D_dB = P_2D_dB - max(P_2D_dB(:));  % Normalizacja do 0 dB
% 
%         figure('Position', [100 100 1000 700]);
%         imagesc(rad2deg(sort(az_unique)), rad2deg(sort(el_unique)), P_2D_dB);
%         colorbar;
%         xlabel('Azimuth [°]');
%         ylabel('Elevation [°]');
%         title('MUSIC-2D Spectrum [dB]');
%         axis xy;
%         caxis([-30 0]);  % Pokazuj tylko top 30 dB
%         hold on;
% 
%         % Wykryte
%         plot(detected_az, detected_el, 'rx', 'MarkerSize', 25, 'LineWidth', 4);
% 
%         % Prawdziwe
%         plot(true_az, true_el, 'go', 'MarkerSize', 20, 'LineWidth', 3);
% 
%         legend('', 'Wykryte', 'Prawdziwe', 'Location', 'best');
%         %grid on;
%     end
% end