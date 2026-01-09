clc; clear; close all;

%% -- config i parametry
config.n_tests = 2;                   % liczba pozycji źródła do przetestowania
config.snr_db = [25, 40, 75];                 % różne poziomy SNR
%config.show_plots_for_test = 3;       % dla której iteracji pokazać wykresy (0 = żadne)

sys.C = 343;
sys.f_min = 300;
sys.f_max = 3300;
sys.f_0 = 3000;
sys.f_cutoff = 4000;
sys.f_s = 16000;
sys.L = 1024;
sys.d = 0.042;
sys.M_row_line = 4;
sys.M = sys.M_row_line * sys.M_row_line;
sys.if_reverbed = 0; % << 0 - false, 1 - true

%% -- mikrofony, macierz i akwizycja
% jak wcześniej, URA 4x4 co 4.2 cm, ICS-41350
mic_element = phased.OmnidirectionalMicrophoneElement("FrequencyRange", [20 10000]);
mic_array = phased.URA("Element", mic_element, "ElementSpacing", sys.d, ...
    "Size", [sys.M_row_line sys.M_row_line]);

mic_positions = getElementPosition(mic_array);
mic_array_centre = mean(mic_positions, 2);

% okno Czebyszewa - poprawia kierunkowość i odpowiedź macierzy
taper = chebwin(sys.M_row_line, 30);
taper2D = taper * taper';
mic_array.Taper = taper2D;

mic_normals = zeros(2, sys.M);
mic_array_conf = phased.ConformalArray('Element', mic_element, ...
    'ElementPosition', mic_positions, 'ElementNormal', mic_normals);
wide_collector = phased.WidebandCollector("Sensor", mic_array_conf, ...
    "Wavefront", "Plane", "PropagationSpeed", sys.C, ...
    "ModulatedInput", false, 'SampleRate', sys.f_s, 'Polarization', 'None');

%% -- grid setup
az_grid = deg2rad(-85:5:85);
el_grid = deg2rad(-45:5:45);
[az_mesh, el_mesh] = meshgrid(az_grid, el_grid);
G_deg_pairs = [az_mesh(:), el_mesh(:)];
G_spherical_vec = az_el_to_spherical_vec(G_deg_pairs(:, 1), G_deg_pairs(:, 2));

% wsparcie dla LUTów
srp_precomp = srp_precompute(mic_positions, G_spherical_vec, ...
    sys.f_min, sys.f_max, sys.f_s, sys.L, sys.C);

%% -- przygotowanie źródeł
rng(41);
room_dim = [4; 4; 3];
test_sources = zeros(3, config.n_tests);
% for i = 1:config.n_tests
%     test_sources(:, i) = [
%         rand() * 2 + 1;
%         rand() * 2 - 1;
%         rand() + 0.5;
%     ];
% end

% -- for n_tests == 2
test_sources(:, 1) = [1; 0; 0.8]; 
test_sources(:, 2) = [2; -0.8; 0.23]; 

%% -- testowy sygnał
% << TODO >> sygnał biały ze spektrm mowy
% t = (0:(1/sys.f_s):((sys.L-1)/sys.f_s))';
% freqs_test = [350, 700, 999, 1111, 3000];
% sig_base = zeros(length(t), 1);
% for f = freqs_test
%     sig_base = sig_base + sin(2 * pi * f * t);
% end
% biały szum z zakresu mowy
t = (0:(1/sys.f_s):((sys.L-1)/sys.f_s))';
sig_base = randn(length(t), 1);
sig_base = sig_base / norm(sig_base);
sig_pow = var(sig_base);

%% -- testy
% przygotowanie
algorithms = {'DAS', 'SRP-PHAT', 'MUSIC', 'MVDR'};
n_algs = length(algorithms);
n_snr = length(config.snr_db);

results = struct();
for a = 1:n_algs
    for s = 1:n_snr
        results(a, s).az_errors = zeros(config.n_tests, 1);
        results(a, s).el_errors = zeros(config.n_tests, 1);
        results(a, s).angular_errors = zeros(config.n_tests, 1);
        results(a, s).times = zeros(config.n_tests, 1);
    end
end

if (sys.if_reverbed == 1) 
    fprintf('DoA tests running - reverb ON\n');
else
    fprintf('DoA tests running - no reverb\n');
end

for snr_idx = 1:n_snr
    snr_db = config.snr_db(snr_idx);
    snr_linear = 10^(snr_db/10);
    noise_var = sig_pow / snr_linear;
    
    fprintf('\n--- SNR = %d dB ---\n', snr_db);
    
    for test_idx = 1:config.n_tests
        fprintf('  test %d/%d: ', test_idx, config.n_tests);
        
        source_pos = test_sources(:, test_idx);
        
        % oczekiwane kąty
        pos_diff = source_pos - mic_array_centre;
        az_true = rad2deg(atan2(pos_diff(2), pos_diff(1)));
        el_true = rad2deg(atan2(pos_diff(3), sqrt(pos_diff(1)^2 + pos_diff(2)^2)));
        
        fprintf('src=[%.2f,%.2f,%.2f], true az=%.1f, el=%.1f\n', ...
            source_pos(1), source_pos(2), source_pos(3), az_true, el_true);
        
        ura_rx = wide_collector(sig_base, [az_true; el_true]);

        if (sys.if_reverbed == 1) 
            % dla braku pogłosu należy zmienić na początku skryptu na: 
            % sys.if_reverbed = 0;

            % -- wysoki pogłos --
            room_resp = acousticRoomResponse(room_dim', source_pos', ...
                mic_positions', "SampleRate", sys.f_s, "SoundSpeed", sys.C, ...
                'MaterialAbsorption', [0.1 0.1 0.1 0.1 0.15 0.2]');

            % -- niski pogłos --
            % room_resp = acousticRoomResponse(room_dim', source_pos', ...
            %     mic_positions', "SampleRate", sys.f_s, "SoundSpeed", sys.C, ...
            %     'MaterialAbsorption', [0.3 0.3 0.3 0.3 0.25 0.4]');

            for m = 1:(sys.M)
                ura_rx(:, m) = filter(room_resp(m, :)', 1, ura_rx(:, m));
            end
        end
        
        noise = sqrt(noise_var) * randn(size(ura_rx));
        ura_rx_noisy = ura_rx + noise;
        
        % LPF pro forma
        ura_rx_processed = lowpass(ura_rx_noisy, sys.f_cutoff, sys.f_s);
        
        % ---- algo tutaj ----
        
        % DaS
        das_pre = das_precomp(mic_positions, ura_rx_processed, ... 
                              G_spherical_vec, G_deg_pairs, sys.C);

        [avg_t, ~] = timing_wrapper(@() asl_das_doa_est(ura_rx_processed, ...
                     G_deg_pairs, sys.f_s, das_pre), 5, 2);
        results(1, snr_idx).times(test_idx) = avg_t;
        
        [das_az, das_el, das_info] = asl_das_doa_est(ura_rx_processed, ...
                     G_deg_pairs, sys.f_s, das_pre);

        results(1, snr_idx).az_errors(test_idx) = das_az - az_true;
        results(1, snr_idx).el_errors(test_idx) = das_el - el_true;
        results(1, snr_idx).angular_errors(test_idx) = ...
            angular_distance(az_true, el_true, das_az, das_el);
        
        %if test_idx == config.show_plots_for_test && snr_idx == 2
            % das_spectrum = zeros(size(G_deg_pairs, 1), 1);
            % das_spectrum(das_info.maxIndex) = das_info.maxPower;
        %end
        
        % SRP-PHAT
        [avg_t, ~] = timing_wrapper(@() asl_srp_phat(ura_rx_processed, ...
            mic_positions, sys.f_min, sys.f_max, sys.f_s, ...
            G_deg_pairs, G_spherical_vec, sys.C, sys.L, srp_precomp), 5, 2);
        results(2, snr_idx).times(test_idx) = avg_t;
        
        [srp_az, srp_el, srp_info] = asl_srp_phat(ura_rx_processed, ...
            mic_positions, sys.f_min, sys.f_max, sys.f_s, ...
            G_deg_pairs, G_spherical_vec, sys.C, sys.L, srp_precomp);

        results(2, snr_idx).az_errors(test_idx) = srp_az - az_true;
        results(2, snr_idx).el_errors(test_idx) = srp_el - el_true;
        results(2, snr_idx).angular_errors(test_idx) = ...
            angular_distance(az_true, el_true, srp_az, srp_el);
        
        % MUSIC(-2D)
        lambda_0 = sys.C / sys.f_0;
        [avg_t, ~] = timing_wrapper(@() asl_music_2d(ura_rx_processed, ...
            lambda_0, mic_positions, 1, G_deg_pairs, G_spherical_vec, ...
            sys.f_s), 5, 2);
        results(3, snr_idx).times(test_idx) = avg_t;
        
        [music_az, music_el, music_info] = asl_music_2d(ura_rx_processed, ...
            lambda_0, mic_positions, 1, G_deg_pairs, G_spherical_vec, ...
            sys.f_s);

        results(3, snr_idx).az_errors(test_idx) = music_az - az_true;
        results(3, snr_idx).el_errors(test_idx) = music_el - el_true;
        results(3, snr_idx).angular_errors(test_idx) = ...
            angular_distance(az_true, el_true, music_az, music_el);
        
        % MVDR
        %lambda_1 = sys.C / 250;
        [avg_t, ~] = timing_wrapper(@() asl_mvdr_doa_est(ura_rx_processed, ...
            mic_positions, G_deg_pairs, G_spherical_vec, sys.f_s, lambda_0), 5, 2);
        results(4, snr_idx).times(test_idx) = avg_t;
        
        [mvdr_az, mvdr_el, mvdr_spec] = asl_mvdr_doa_est(ura_rx_processed, ...
            mic_positions, G_deg_pairs, G_spherical_vec, sys.f_s, lambda_0);

        results(4, snr_idx).az_errors(test_idx) = mvdr_az - az_true;
        results(4, snr_idx).el_errors(test_idx) = mvdr_el - el_true;
        results(4, snr_idx).angular_errors(test_idx) = ...
            angular_distance(az_true, el_true, mvdr_az, mvdr_el);
        
        fprintf('    DAS:      az=%.1f /err=%.1f, El=%.1f /err=%.1f\n', ...
            das_az, das_az-az_true, das_el, das_el-el_true);
        fprintf('    SRP-PHAT: az=%.1f /err=%.1f, El=%.1f /err=%.1f\n', ...
            srp_az, srp_az-az_true, srp_el, srp_el-el_true);
        fprintf('    MUSIC:    az=%.1f /err=%.1f, El=%.1f /err=%.1f\n', ...
            music_az, music_az-az_true, music_el, music_el-el_true);
        fprintf('    MVDR:     az=%.1f /err=%.1f, El=%.1f /err=%.1f\n', ...
            mvdr_az, mvdr_az-az_true, mvdr_el, mvdr_el-el_true);
        
        % wykresy
        %if test_idx == config.show_plots_for_test && snr_idx == 2

            plot_data.az_grid = rad2deg(az_grid);
            plot_data.el_grid = rad2deg(el_grid);
            plot_data.az_true = az_true;
            plot_data.el_true = el_true;
            plot_data.estimates = struct(...
                'das', struct('az', das_az, 'el', das_el, 'spectrum', das_info.powerSpectrum), ...
                'srp', struct('az', srp_az, 'el', srp_el, 'spectrum', srp_info.srp_map), ...
                'music', struct('az', music_az, 'el', music_el, 'spectrum', music_info.P_MUSIC), ...
                'mvdr', struct('az', mvdr_az, 'el', mvdr_el, 'spectrum', mvdr_spec));
            
            plot_doa_results(plot_data, G_deg_pairs, test_idx, snr_db);
        %end
    end
end

%% -- tutaj statystyki i wykresy końcowe
fprintf('\n========================================\n');

plot_comparative_results(results, algorithms, config.snr_db);

for snr_idx = 1:n_snr
    fprintf('\n=== SNR = %d dB ===\n', config.snr_db(snr_idx));
    fprintf('%-12s | %8s | %8s | %10s | %10s\n', ...
        'Algorytm', 'Az RMSE', 'El RMSE', 'Ang RMSE', 'Czas[ms]');
    fprintf('-----------------------------------------------------------\n');
    
    for a = 1:n_algs
        az_rmse = sqrt(mean(results(a, snr_idx).az_errors.^2));
        el_rmse = sqrt(mean(results(a, snr_idx).el_errors.^2));
        ang_rmse = sqrt(mean(results(a, snr_idx).angular_errors.^2));
        time_avg = mean(results(a, snr_idx).times) * 1000;
        
        fprintf('%-12s | %7.2f | %7.2f | %9.2f | %8.2f\n', ...
            algorithms{a}, az_rmse, el_rmse, ang_rmse, time_avg);
    end
end

%% ==== funkcje pomocnicze ====

function ang_dist = angular_distance(az1, el1, az2, el2)
    az1_rad = deg2rad(az1);
    el1_rad = deg2rad(el1);
    az2_rad = deg2rad(az2);
    el2_rad = deg2rad(el2);
    
    v1 = [cos(el1_rad)*cos(az1_rad); cos(el1_rad)*sin(az1_rad); sin(el1_rad)];
    v2 = [cos(el2_rad)*cos(az2_rad); cos(el2_rad)*sin(az2_rad); sin(el2_rad)];
    
    cos_angle = dot(v1, v2);
    cos_angle = max(-1, min(1, cos_angle));
    ang_dist = rad2deg(acos(cos_angle));
end

function plot_doa_results(data, G_deg_pairs, test_num, snr_db)
    figure('Position', [100, 100, 1600, 900]);
    sgtitle(sprintf('Mapy mocy algorytmów DoA (Test #%d, SNR=%d dB)', test_num, snr_db), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    algorithms = {'DAS', 'SRP-PHAT', 'MUSIC', 'MVDR'};
    spectra = {data.estimates.das.spectrum, data.estimates.srp.spectrum, ...
                data.estimates.music.spectrum, data.estimates.mvdr.spectrum};
    estimates = {data.estimates.das, data.estimates.srp, ...
                 data.estimates.music, data.estimates.mvdr};
    
    for i = 1:4
        subplot(2, 2, i);
        
        if ~isempty(spectra{i})
            spectrum_2d = reshape(spectra{i}, length(data.el_grid), length(data.az_grid));
            % --- przesunięcie, żeby wszystkie wartości były ≥ 0
            spectrum_shifted = spectrum_2d - min(spectrum_2d(:));
            % --- normalizacja do [0,1]
            spectrum_shifted = spectrum_shifted / max(spectrum_shifted(:));
            % --- log dB
            spectrum_db = 10*log10(spectrum_shifted + eps); % eps, żeby uniknąć log(0)
            
            imagesc(data.az_grid, data.el_grid, spectrum_db);
            colorbar;
            caxis([-10 0]);
        end
        
        hold on;

        plot(data.az_true, data.el_true, 'go', 'MarkerSize', 12, ...
            'LineWidth', 3, 'DisplayName', 'Prawdziwa');

        plot(estimates{i}.az, estimates{i}.el, 'rx', 'MarkerSize', 15, ...
            'LineWidth', 3, 'DisplayName', 'Estymacja');

        hold off;
        
        xlabel('Azymut [°]');
        ylabel('Elewacja [°]');
        title(algorithms{i});
        legend('Location', 'best');
        axis xy;
        grid on;
    end
end

function plot_comparative_results(results, algorithms, snr_levels)
    figure('Position', [100, 100, 1400, 500]);
    
    n_algs = length(algorithms);
    n_snr = length(snr_levels);
    
    ang_rmse = zeros(n_algs, n_snr);
    for a = 1:n_algs
        for s = 1:n_snr
            ang_rmse(a, s) = sqrt(mean(results(a, s).angular_errors.^2));
        end
    end
    
    subplot(1, 3, 1);
    colors = lines(n_algs);
    for a = 1:n_algs
        plot(snr_levels, ang_rmse(a, :), '-o', 'LineWidth', 2, ...
            'Color', colors(a, :), 'MarkerSize', 8, 'DisplayName', algorithms{a});
        hold on;
    end
    hold off;
    xlabel('SNR [dB]');
    ylabel('Średni błąd kątowy [°]');
    title('Dokładność vs SNR');
    legend('Location', 'best');
    grid on;
    
    subplot(1, 3, 2);
    mid_snr_idx = ceil(n_snr / 2);
    
    hold on;
    box_width = 0.6;
    for a = 1:n_algs
        errors = results(a, mid_snr_idx).angular_errors;
        
        q1 = prctile(errors, 25);
        q2 = median(errors);
        q3 = prctile(errors, 75);
        iqr = q3 - q1;
        lower_whisker = max(min(errors), q1 - 1.5*iqr);
        upper_whisker = min(max(errors), q3 + 1.5*iqr);
        outliers = errors(errors < lower_whisker | errors > upper_whisker);
        
        x = a;
        rectangle('Position', [x-box_width/2, q1, box_width, q3-q1], ...
            'FaceColor', colors(a,:), 'EdgeColor', 'k', 'LineWidth', 1.5);
        
        plot([x-box_width/2, x+box_width/2], [q2, q2], 'k-', 'LineWidth', 2);
        
        plot([x, x], [q3, upper_whisker], 'k-', 'LineWidth', 1.5);
        plot([x, x], [lower_whisker, q1], 'k-', 'LineWidth', 1.5);
        plot([x-box_width/4, x+box_width/4], [upper_whisker, upper_whisker], 'k-', 'LineWidth', 1.5);
        plot([x-box_width/4, x+box_width/4], [lower_whisker, lower_whisker], 'k-', 'LineWidth', 1.5);
        
        if ~isempty(outliers)
            plot(repmat(x, size(outliers)), outliers, 'r+', 'MarkerSize', 6);
        end
    end
    hold off;
    
    xlim([0.5, n_algs+0.5]);
    xticks(1:n_algs);
    xticklabels(algorithms);
    ylabel('Błąd kątowy [°]');
    title(sprintf('Rozkład błędów (SNR = %d dB)', snr_levels(mid_snr_idx)));
    grid on;
    
    % błędy
    subplot(1, 3, 3);
    hold on;
    for a = 1:n_algs
        errors = results(a, mid_snr_idx).angular_errors;
        histogram(errors, 'FaceColor', colors(a,:), 'FaceAlpha', 0.5, ...
            'DisplayName', algorithms{a}, 'BinWidth', 0.5);
    end
    hold off;
    xlabel('Błąd kątowy [°]');
    ylabel('Liczba wystąpień');
    title(sprintf('Histogram błędów (SNR = %d dB)', snr_levels(mid_snr_idx)));
    legend('Location', 'best');
    grid on;
end