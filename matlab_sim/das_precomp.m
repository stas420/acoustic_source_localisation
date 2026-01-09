function [das_taus] = das_precomp(micsPositions, rxSignal, ...
                                  gridDirVectors, gridAnglePairs, C)
    G = size(gridAnglePairs, 1);
    [~, M] = size(rxSignal);

    das_taus = zeros(G, M);
    mic_centre = mean(micsPositions, 2);

    %% -- precalc of TDoAs
    for m = 1:M
        d = (micsPositions(:, m) - mic_centre);  % [3 x 1]
        
        % gridDirVectors: [3 x G]
        % d: [3 x 1]
        % [G x 3] * [3 x 1] = [G x 1]
        das_taus(:, m) = (gridDirVectors' * d) / C;
    end
end