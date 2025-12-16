%% Framework do testowania algorytmów DoA
% Przeprowadza serię testów i zbiera statystyki błędów
clc; clear; close all;

%% ===== KONFIGURACJA TESTÓW =====
config.n_tests = 10;                    % liczba pozycji źródła do przetestowania
config.snr_db = [10, 20, 30];          % różne poziomy SNR do przetestowania
config.show_plots_for_test = 5;        % dla której iteracji pokazać wykresy (0 = nie pokazuj)

%% ===== PARAMETRY SYSTEMU (stałe) =====
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

%% ===== SETUP MIKROFONU =====
mic_element = phased.OmnidirectionalMicrophoneElement("FrequencyRange", [20 10000]);
mic_array = phased.URA("Element", mic_element, "ElementSpacing", sys.d, ...
    "Size", [sys.M_row_line sys.M_row_line]);
mic_positions = getElementPosition(mic_array);
mic_array_centre = mean(mic_positions, 2);

% Okno Chebysheva
taper = chebwin(sys.M_row_line, 20);
taper2D = taper * taper';
mic_array.Taper = taper2D;

% WidebandCollector
mic_normals = zeros(2, sys.M);
mic_array_conf = phased.ConformalArray('Element', mic_element, ...
    'ElementPosition', mic_positions, 'ElementNormal', mic_normals);
wide_collector = phased.WidebandCollector("Sensor", mic_array_conf, ...
    "Wavefront", "Plane", "PropagationSpeed", sys.C, ...
    "ModulatedInput", false, 'SampleRate', sys.f_s, 'Polarization', 'None');

%% ===== SIATKA KĄTÓW =====
az_grid = deg2rad(-85:2.5:85);  % gęstsza siatka dla lepszej precyzji
el_grid = deg2rad(-45:2.5:45);
[az_mesh, el_mesh] = meshgrid(az_grid, el_grid);
G_deg_pairs = [az_mesh(:), el_mesh(:)];
G_spherical_vec = az_el_to_spherical_vec(G_deg_pairs(:, 1), G_deg_pairs(:, 2));

%% ===== PRE-COMPUTE DLA SRP-PHAT =====
srp_precomp = srp_precompute(mic_positions, G_spherical_vec, ...
    sys.f_min, sys.f_max, sys.f_s, sys.L, sys.C);

%% ===== GENEROWANIE POZYCJI ŹRÓDEŁ =====
% Losowe pozycje w rozsądnym zakresie
rng(42); % dla powtarzalności
room_dim = [4; 4; 3];
test_sources = zeros(3, config.n_tests);
for i = 1:config.n_tests
    % Losowe pozycje w pokoju, unikając kątów skrajnych
    test_sources(:, i) = [
        rand() * 3 - 1.5;  % x: -1.5 do 1.5
        rand() * 2 - 1;    % y: -1 do 1
        rand() * 1.5 + 0.3 % z: 0.3 do 1.8
    ];
end

%% ===== PRZYGOTOWANIE SYGNAŁU TESTOWEGO =====
t = (0:(1/sys.f_s):((sys.L-1)/sys.f_s))';
freqs_test = [350, 700, 999, 1111, 3000];
sig_base = zeros(length(t), 1);
for f = freqs_test
    sig_base = sig_base + sin(2 * pi * f * t);
end
sig_pow = var(sig_base);

%% ===== INICJALIZACJA WYNIKÓW =====
algorithms = {'DAS', 'SRP-PHAT', 'MUSIC', 'MVDR'};
n_algs = length(algorithms);
n_snr = length(config.snr_db);

% Struktura wyników: results(alg_idx, snr_idx)
results = struct();
for a = 1:n_algs
    for s = 1:n_snr
        results(a, s).az_errors = zeros(config.n_tests, 1);
        results(a, s).el_errors = zeros(config.n_tests, 1);
        results(a, s).angular_errors = zeros(config.n_tests, 1);
        results(a, s).times = zeros(config.n_tests, 1);
    end
end

%% ===== PĘTLA TESTOWA =====
fprintf('========================================\n');
fprintf('ROZPOCZYNAM TESTY DoA\n');
fprintf('Liczba testów: %d\n', config.n_tests);
fprintf('Poziomy SNR: %s dB\n', num2str(config.snr_db));
fprintf('========================================\n\n');

for snr_idx = 1:n_snr
    snr_db = config.snr_db(snr_idx);
    snr_linear = 10^(snr_db/10);
    noise_var = sig_pow / snr_linear;
    
    fprintf('\n--- SNR = %d dB ---\n', snr_db);
    
    for test_idx = 1:config.n_tests
        fprintf('  Test %d/%d: ', test_idx, config.n_tests);
        
        % Pozycja źródła
        source_pos = test_sources(:, test_idx);
        
        % Prawdziwe kąty
        pos_diff = source_pos - mic_array_centre;
        az_true = rad2deg(atan2(pos_diff(2), pos_diff(1)));
        el_true = rad2deg(atan2(pos_diff(3), sqrt(pos_diff(1)^2 + pos_diff(2)^2)));
        
        fprintf('Src=[%.2f,%.2f,%.2f], True Az=%.1f°, El=%.1f°\n', ...
            source_pos(1), source_pos(2), source_pos(3), az_true, el_true);
        
        % Akwizycja sygnału
        ura_rx = wide_collector(sig_base, [az_true; el_true]);
        
        % Dodanie szumu
        noise = sqrt(noise_var) * randn(size(ura_rx));
        ura_rx_noisy = ura_rx + noise;
        
        % LPF
        ura_rx_processed = lowpass(ura_rx_noisy, sys.f_cutoff, sys.f_s);
        
        %% === ALGORYTMY ===
        
        % 1. DAS
        [avg_t, ~] = accurate_timing_wrapper(@() asl_das_doa_est(ura_rx_processed, ...
            mic_positions, G_deg_pairs, G_spherical_vec, sys.f_s), 5, 2);
        results(1, snr_idx).times(test_idx) = avg_t;
        
        [das_az, das_el, das_info] = asl_das_doa_est(ura_rx_processed, ...
            mic_positions, G_deg_pairs, G_spherical_vec, sys.f_s);
        % Normalizacja kątów
        [das_az, das_el] = normalize_doa_angles(das_az, das_el, az_true, el_true);
        results(1, snr_idx).az_errors(test_idx) = das_az - az_true;
        results(1, snr_idx).el_errors(test_idx) = das_el - el_true;
        results(1, snr_idx).angular_errors(test_idx) = ...
            angular_distance(az_true, el_true, das_az, das_el);
        if test_idx == config.show_plots_for_test && snr_idx == 2
            das_spectrum = zeros(size(G_deg_pairs, 1), 1);
            das_spectrum(das_info.maxIndex) = das_info.maxPower;
        end
        
        % 2. SRP-PHAT
        [avg_t, ~] = accurate_timing_wrapper(@() asl_srp_phat(ura_rx_processed, ...
            mic_positions, sys.f_min, sys.f_max, sys.f_s, ...
            G_deg_pairs, G_spherical_vec, sys.C, sys.L, srp_precomp), 5, 2);
        results(2, snr_idx).times(test_idx) = avg_t;
        
        [srp_az, srp_el, srp_info] = asl_srp_phat(ura_rx_processed, ...
            mic_positions, sys.f_min, sys.f_max, sys.f_s, ...
            G_deg_pairs, G_spherical_vec, sys.C, sys.L, srp_precomp);
        % Normalizacja kątów
        [srp_az, srp_el] = normalize_doa_angles(srp_az, srp_el, az_true, el_true);
        results(2, snr_idx).az_errors(test_idx) = srp_az - az_true;
        results(2, snr_idx).el_errors(test_idx) = srp_el - el_true;
        results(2, snr_idx).angular_errors(test_idx) = ...
            angular_distance(az_true, el_true, srp_az, srp_el);
        
        % 3. MUSIC
        lambda_0 = sys.C / sys.f_0;
        [avg_t, ~] = accurate_timing_wrapper(@() asl_music_2d(ura_rx_processed, ...
            lambda_0, sys.d, sys.M_row_line, 1, G_deg_pairs), 5, 2);
        results(3, snr_idx).times(test_idx) = avg_t;
        
        [music_az, music_el, music_info] = asl_music_2d(ura_rx_processed, ...
            lambda_0, sys.d, sys.M_row_line, 1, G_deg_pairs);
        % Normalizacja kątów
        [music_az, music_el] = normalize_doa_angles(music_az, music_el, az_true, el_true);
        results(3, snr_idx).az_errors(test_idx) = music_az - az_true;
        results(3, snr_idx).el_errors(test_idx) = music_el - el_true;
        results(3, snr_idx).angular_errors(test_idx) = ...
            angular_distance(az_true, el_true, music_az, music_el);
        
        % 4. MVDR
        [avg_t, ~] = accurate_timing_wrapper(@() asl_mvdr_doa_est(ura_rx_processed, ...
            mic_positions, G_deg_pairs, G_spherical_vec, sys.f_s, lambda_0), 5, 2);
        results(4, snr_idx).times(test_idx) = avg_t;
        
        [mvdr_az, mvdr_el, mvdr_spec] = asl_mvdr_doa_est(ura_rx_processed, ...
            mic_positions, G_deg_pairs, G_spherical_vec, sys.f_s, lambda_0);
        results(4, snr_idx).az_errors(test_idx) = mvdr_az - az_true;
        results(4, snr_idx).el_errors(test_idx) = mvdr_el - el_true;
        results(4, snr_idx).angular_errors(test_idx) = ...
            angular_distance(az_true, el_true, mvdr_az, mvdr_el);
        
        % Wyświetl wyniki dla tego testu
        fprintf('    DAS:      Az=%.1f° (err=%.1f°), El=%.1f° (err=%.1f°)\n', ...
            das_az, das_az-az_true, das_el, das_el-el_true);
        fprintf('    SRP-PHAT: Az=%.1f° (err=%.1f°), El=%.1f° (err=%.1f°)\n', ...
            srp_az, srp_az-az_true, srp_el, srp_el-el_true);
        fprintf('    MUSIC:    Az=%.1f° (err=%.1f°), El=%.1f° (err=%.1f°)\n', ...
            music_az, music_az-az_true, music_el, music_el-el_true);
        fprintf('    MVDR:     Az=%.1f° (err=%.1f°), El=%.1f° (err=%.1f°)\n', ...
            mvdr_az, mvdr_az-az_true, mvdr_el, mvdr_el-el_true);
        
        %% === WYKRESY DLA WYBRANEGO TESTU ===
        if test_idx == config.show_plots_for_test && snr_idx == 2
            % Przygotuj dane do plotowania
            plot_data.az_grid = rad2deg(az_grid);
            plot_data.el_grid = rad2deg(el_grid);
            plot_data.az_true = az_true;
            plot_data.el_true = el_true;
            plot_data.estimates = struct(...
                'das', struct('az', das_az, 'el', das_el, 'spectrum', das_info.maxPower), ...
                'srp', struct('az', srp_az, 'el', srp_el, 'spectrum', []), ...
                'music', struct('az', music_az, 'el', music_el, 'spectrum', music_info.P_MUSIC), ...
                'mvdr', struct('az', mvdr_az, 'el', mvdr_el, 'spectrum', mvdr_spec));
            
            % Wywołaj funkcję plotującą
            plot_doa_results(plot_data, G_deg_pairs, test_idx, snr_db);
        end
    end
end

%% ===== STATYSTYKI I WYNIKI =====
fprintf('\n========================================\n');
fprintf('PODSUMOWANIE WYNIKÓW\n');
fprintf('========================================\n\n');

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
        
        fprintf('%-12s | %7.2f° | %7.2f° | %9.2f° | %8.2f\n', ...
            algorithms{a}, az_rmse, el_rmse, ang_rmse, time_avg);
    end
end

%% ===== WYKRESY PORÓWNAWCZE =====
plot_comparative_results(results, algorithms, config.snr_db);

%% ===== FUNKCJE POMOCNICZE =====
function ang_dist = angular_distance(az1, el1, az2, el2)
    % Oblicza odległość kątową między dwoma punktami na sferze
    az1_rad = deg2rad(az1);
    el1_rad = deg2rad(el1);
    az2_rad = deg2rad(az2);
    el2_rad = deg2rad(el2);
    
    % Konwersja do wektorów jednostkowych
    v1 = [cos(el1_rad)*cos(az1_rad); cos(el1_rad)*sin(az1_rad); sin(el1_rad)];
    v2 = [cos(el2_rad)*cos(az2_rad); cos(el2_rad)*sin(az2_rad); sin(el2_rad)];
    
    % Odległość kątowa przez iloczyn skalarny
    cos_angle = dot(v1, v2);
    cos_angle = max(-1, min(1, cos_angle)); % clamp dla stabilności numerycznej
    ang_dist = rad2deg(acos(cos_angle));
end

function plot_doa_results(data, G_deg_pairs, test_num, snr_db)
    figure('Position', [100, 100, 1600, 900]);
    sgtitle(sprintf('Mapy mocy algorytmów DoA (Test #%d, SNR=%d dB)', test_num, snr_db), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    algorithms = {'DAS', 'SRP-PHAT', 'MUSIC', 'MVDR'};
    spectra = {[], [], data.estimates.music.spectrum, data.estimates.mvdr.spectrum};
    estimates = {data.estimates.das, data.estimates.srp, ...
                 data.estimates.music, data.estimates.mvdr};
    
    for i = 1:4
        subplot(2, 2, i);
        
        if ~isempty(spectra{i})
            % Reshape do 2D
            spectrum_2d = reshape(spectra{i}, length(data.el_grid), length(data.az_grid));
            spectrum_db = 10*log10(spectrum_2d / max(spectrum_2d(:)));
            
            imagesc(data.az_grid, data.el_grid, spectrum_db);
            colorbar;
            caxis([-30 0]);
        else
            % Dla DAS i SRP-PHAT nie mamy pełnego spektrum
            text(0.5, 0.5, 'Spektrum niedostępne', ...
                'HorizontalAlignment', 'center', 'Units', 'normalized');
        end
        
        hold on;
        % Prawdziwa pozycja (zielone kółko)
        plot(data.az_true, data.el_true, 'go', 'MarkerSize', 12, ...
            'LineWidth', 3, 'DisplayName', 'Prawdziwa');
        % Estymacja (czerwony krzyżyk)
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
    % Wykres porównawczy RMSE vs SNR
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
    
    % Wykres pudełkowy własnej roboty dla błędów przy średnim SNR
    subplot(1, 3, 2);
    mid_snr_idx = ceil(n_snr / 2);
    
    % Ręczna implementacja boxplot
    hold on;
    box_width = 0.6;
    for a = 1:n_algs
        errors = results(a, mid_snr_idx).angular_errors;
        
        % Statystyki
        q1 = prctile(errors, 25);
        q2 = median(errors);
        q3 = prctile(errors, 75);
        iqr = q3 - q1;
        lower_whisker = max(min(errors), q1 - 1.5*iqr);
        upper_whisker = min(max(errors), q3 + 1.5*iqr);
        outliers = errors(errors < lower_whisker | errors > upper_whisker);
        
        % Rysuj pudełko
        x = a;
        rectangle('Position', [x-box_width/2, q1, box_width, q3-q1], ...
            'FaceColor', colors(a,:), 'EdgeColor', 'k', 'LineWidth', 1.5);
        
        % Mediana
        plot([x-box_width/2, x+box_width/2], [q2, q2], 'k-', 'LineWidth', 2);
        
        % Wąsy
        plot([x, x], [q3, upper_whisker], 'k-', 'LineWidth', 1.5);
        plot([x, x], [lower_whisker, q1], 'k-', 'LineWidth', 1.5);
        plot([x-box_width/4, x+box_width/4], [upper_whisker, upper_whisker], 'k-', 'LineWidth', 1.5);
        plot([x-box_width/4, x+box_width/4], [lower_whisker, lower_whisker], 'k-', 'LineWidth', 1.5);
        
        % Outliery
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
    
    % Dodatkowy wykres: histogram błędów
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