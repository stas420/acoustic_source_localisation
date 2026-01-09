function [avg_time, std_time] = timing_wrapper(func_handle, n_runs, warmup_runs)
    % this function is for proper timing checks:
    % - removes compilation oerhead impact
    % - makes multiple runs and averages
    % - prevents Matlab from reusing pre-calculated data (more like in real-time)
    
    if nargin < 2, n_runs = 10; end
    if nargin < 3, warmup_runs = 3; end
    
    % run a few time to have JIT prepared
    for i = 1:warmup_runs
        func_handle();
    end
    
    times = zeros(n_runs, 1);
    
    for i = 1:n_runs
        % GC run to have clear stage
        java.lang.System.gc();
        pause(0.01); % short break, GC takes time
        
        t_start = tic;
        func_handle();
        times(i) = toc(t_start);
    end
    
    % remove outliers
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