function [avg_time, std_time] = accurate_timing_wrapper(func_handle, n_runs, warmup_runs)
    % DOKŁADNY POMIAR CZASU WYKONANIA FUNKCJI
    % Eliminuje wpływ JIT compilation, cache, i innych artefaktów
    %
    % Inputs:
    %   func_handle - funkcja do zmierzenia, np: @() my_function(args)
    %   n_runs - liczba powtórzeń do uśrednienia (domyślnie 10)
    %   warmup_runs - liczba przebiegów rozgrzewkowych (domyślnie 3)
    %
    % Outputs:
    %   avg_time - średni czas [sekundy]
    %   std_time - odchylenie standardowe [sekundy]
    
    if nargin < 2, n_runs = 10; end
    if nargin < 3, warmup_runs = 3; end
    
    % === WARMUP: uruchom kilka razy żeby JIT się skompilował ===
    for i = 1:warmup_runs
        func_handle();
    end
    
    % === WŁAŚCIWE POMIARY ===
    times = zeros(n_runs, 1);
    
    for i = 1:n_runs
        % Wymuś garbage collection przed pomiarem
        java.lang.System.gc();
        pause(0.01); % krótka pauza żeby GC się zakończył
        
        % Pomiar
        t_start = tic;
        func_handle();
        times(i) = toc(t_start);
    end
    
    % Usuń outliery (np. z powodu przełączenia kontekstu systemu)
    % Metoda: usuń najwyższy i najniższy wynik jeśli n_runs >= 5
    if n_runs >= 5
        times_sorted = sort(times);
        times_trimmed = times_sorted(2:end-1);
        avg_time = mean(times_trimmed);
        std_time = std(times_trimmed);
    else
        avg_time = mean(times);
        std_time = std(times);
    end
end