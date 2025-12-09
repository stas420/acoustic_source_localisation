%% ---- najwazniejsze zmienne ----
clc; clear;
C = 343;            % prędkość dźwięku
f_min = 300;        % minimalna częstotliwość, która nas interesuje
f_max = 3300;       % maksymalna częstotliwość, która nas interesuje
f_0 = 3000;         % częstotliwość zainteresowania, np. centralna dla narrowband
f_cutoff = 4000;    % częstotliwość graniczna dla LPF
f_s = 16000;        % częstotliwość próbkowania
L = 1024;           % długość ramki danych per mikrofon
d = 0.042;          % odległości w pionie i poziomie między mikrofonami
M_row_line = 4;     % ilość mikrofonów w kolumnie i wierszu
M = M_row_line * M_row_line; % łączna ilość mikrofonów w macierzy

%% ---- macierz, pokój i ich geometrie ----
% wzór mikrofonu: Adafruit PDM MEMS Omni mic 3492 (Botland) 
% -- ważne: tu nie robię filtrowania PDM->PCM, zakładam, że dostaję od razu
% --        gotowy sygnał dźwiękowy do przetwarzania
mic_element = phased.OmnidirectionalMicrophoneElement("FrequencyRange", ...
    [20 10000]);

% ta macierz jest "zawarta" w płaszczyznie yz, a "zwrócona przodem" zgodnie
% z osią x
mic_array = phased.URA("Element", mic_element, "ElementSpacing", d, ...
    "Size", [M_row_line M_row_line]);
mic_positions_relative = getElementPosition(mic_array);
%mic_array_change = [1; 1; 1];
mic_positions_absolute = mic_positions_relative; % + mic_array_change;
mic_array_centre = mean(mic_positions_absolute, 2);

% punkt odniesienia, względem którego liczone będą [az;el] znajduje się
% w środku macierzy, a on znajduje się w punkcie 0,0,0 (zgodnie z
% dokumentacją phased.URA)
room_dim = [4; 4; 3]; % przykładowe wymiary pokoju
source_pos = [1.9; -0.3; 0.9]; % przykładowa pozycja źródła

mic_normals = zeros(2,M);
mic_normals(1,:) = 0;
mic_normals(2,:) = 0;
mic_array_conf = phased.ConformalArray('Element', mic_element, ...
                   'ElementPosition', mic_positions_absolute, ...
                   'ElementNormal', mic_normals);
mic_array_conf_pos_check = getElementPosition(mic_array_conf);

% mic_positions_absolute(1,:) = mic_positions_relative(1,:) + mic_array_change(1);
% mic_positions_absolute(2,:) = mic_positions_relative(2,:) + mic_array_change(2);
% mic_positions_absolute(3,:) = mic_positions_relative(3,:) + mic_array_change(3);
writematrix(mic_positions_absolute, 'mics.txt'); % export do testów implementacji w C

% okno Chebysheva, by poprawić tłumienie po bokach - bez tego jest trochę słabo
% -- szersza główna wiązka, ale nam to nie przeszkadza
taper = chebwin(M_row_line, 20);   % 20 dB tłumienia na linii 4 mikrofonów
taper2D = taper * taper'; % przejście do 2D
mic_array.Taper = taper2D;% nakładam okno

% "echo" pokoju w formie odpowiedzi impulsowej, czyli de facto filtra FIR
room_resp = acousticRoomResponse(room_dim', source_pos', ...
    mic_positions_absolute', "SampleRate", f_s, "SoundSpeed", C);

%% ---- przygotowanie sygnału i jego obróbki ----
% kąty padania, na razie uproszczone obliczenia
pos_diff = source_pos - mic_array_centre;
az = rad2deg(atan2(pos_diff(2), pos_diff(1)));
el = rad2deg(atan2(pos_diff(3), sqrt(pos_diff(1)^2 + pos_diff(2)^2)));

% WidebandCollector symuluje przesunięcia fazowe na każdym z elementów
% macierzy i nie obcina przy tym pasma
wide_collector = phased.WidebandCollector("Sensor", mic_array_conf, ...
    "Wavefront", "Plane", "PropagationSpeed", C, ...
    "ModulatedInput", false, 'SampleRate', f_s, 'Polarization', 'None');

% przykładowy sygnał, tj. sinusy z zakresu mowy
t = (0:(1/f_s):((L-1)/f_s))';
freqs = [350, 700, 999, 1111, 3000];
sig = zeros(length(t), 1);

for f = freqs
    sig = sig + sin(2 * pi * f * t);
end

% przygotowanie drobnego szumu do dodania później
sig_pow = var(sig);
snr_target_db = 55;
snr_target = 10^(snr_target_db/10);
noise_var = sig_pow/snr_target;
noise = sqrt(noise_var) * randn(size(sig));

% zobrazowanie pozycji mikrofonów względem źródła dźwięku
% figure;
% scatter3([mic_positions_absolute(1,:) source_pos(1)], ...
%          [mic_positions_absolute(2,:) source_pos(2)], ...
%          [mic_positions_absolute(3,:) source_pos(3)]);

%% ---- właściwa akwizycja sygnału i obróbka
% łapiemy sygnał nadchodzący z odpowiedniego kąta
ura_rx = wide_collector(sig, [az; el]);

% przygotowanie do obróbki: nałożenie "echa" i oszumienie
ura_rx_rev = zeros(L, M);
ura_rx_rev_noise = zeros(L, M);

for m = 1:M
    % filter, czyli nakładamy odpowiedź pokoju jako FIR (stąd mianownik 1)
    ura_rx_rev(:, m) = ura_rx(:, m); %filter(room_resp(m, :)', 1, ura_rx(:, m)); << TODO, it destroys everything
    ura_rx_rev_noise(:, m) = ura_rx_rev(:, m) + noise; % << check if it also destroys everything....
end

% LPF pro forma
ura_rx_rev_noise_low = lowpass(ura_rx_rev_noise, f_cutoff, f_s);
writematrix(ura_rx, 'ura_rx.txt'); % export do testów implementacji w C

% wykresy dla porównania
% f_axis = (f_s/L)*(0:(L-1));
% figure;
% subplot(2, 2, 1);
% plot(sig);
% title('sig');
% subplot(2, 2, 2);
% plot(ura_rx_rev_noise_low(:, 1));
% title('ura rx rev noise low');
% subplot(2, 2, 3);
% plot(f_axis, abs(fft(sig)));
% title('fft(sig)');
% subplot(2, 2, 4);
% plot(f_axis, abs(fft(ura_rx_rev_noise_low(:, 1))));
% title('fft(ura rx rev noise low)');

%% ---- algorytmy ASL/SSL ----
% przygotowanie "wspólnych" kwestii
% wspólna siatka kierunków do przeszukiwania
az_grid = deg2rad(-85:5:85);
el_grid = deg2rad(-45:5:45);
[az_mesh, el_mesh] = meshgrid(az_grid, el_grid);
G_deg_pairs = [az_mesh(:), el_mesh(:)];
G_spherical_vec = az_el_to_spherical_vec(G_deg_pairs(:, 1), G_deg_pairs(:, 2));

%% -- delay-and-sum --
[das_az, das_el, das_info] = asl_das_doa_est(ura_rx_rev_noise_low, ...
                             mic_positions_absolute, G_deg_pairs, ...
                             G_spherical_vec, f_s);

%% -- SRP-PHAT -- << done, working
[srp_phat_az, srp_phat_el, srp_phat_info] = asl_srp_phat(ura_rx_rev_noise_low, ... 
    mic_positions_absolute, f_min, f_max, f_s, G_deg_pairs, G_spherical_vec, C, L);

%% -- MUSIC -- << done, almost working - angles conv to be checked
lambda_0 = C/f_0;
[music_az, music_el, music_info] = asl_music_2d(ura_rx_rev_noise_low, ...
    lambda_0, d, M_row_line, 1, G_deg_pairs);

%% -- MVDR --
[mvdr_az, mvdr_el, mvdr_spec] = asl_mvdr_doa_est(ura_rx_rev_noise_low, ...
    mic_positions_absolute, G_deg_pairs, G_spherical_vec, f_s, lambda_0);