function g_phat = asl_gcc_phat(X_l, X_m)
    % GCC-PHAT (Generalized Cross-Correlation with Phase Transform)
    %
    % inputs:
    %   X_l - DFT of signal x received on mic l [n_freqs x 1]
    %   X_m - DFT of signal x received on mic m [n_freqs x 1]
    %
    % output:
    %   g_phat - GCC-PHAT func for each freq [n_freqs x 1]
    %            where each complex value represents phase shift
    %
    % reference equation:
    %   g_phat(f) = [X_l(f) · X_m*(f)] / |X_l(f) · X_m*(f)|
    %
    % PHAT removes amp info and leaves only phase(-shifts) info,
    % as the amp is irrelevant and we only care about the phase
    
    % check if the vectors are [n x 1]
    X_l = X_l(:);
    X_m = X_m(:);
    
    cross_spectrum = X_l .* conj(X_m);
    % gonna be adding epsilon for numeric reasons -- [ref. needed!]
    epsilon = 1e-10;    

    % g_phat is a complex-valued vector (function in freq domain)
    % where each value has amp of 1 and phase equal to phase diff between lm
    %
    %   g_phat(f) = e^(j * delta_phi(f))
    %   where delta phi(f) = phi_l(f) - phi_m(f)
    g_phat = cross_spectrum ./ (abs(cross_spectrum) + epsilon);
end