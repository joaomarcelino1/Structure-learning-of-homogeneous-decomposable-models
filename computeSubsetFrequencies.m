function result = computeSubsetFrequencies(binaryVectors)
    % Input: Binary dataset
    % Output: result - A cell array where each cell corresponds to a subset of coordinates
    %         and contains frequency information for that subset
    [N, n] = size(binaryVectors);
    
    result = cell(2^n, 1);
    
    for i = 0:(2^n - 1)
        % Convert i to binary to represent the subset
        subset = dec2bin(i, n) == '1';
        
        if sum(subset) == 0
            % Empty set - no need to compute frequencies
            result{i+1} = struct('subset', [], 'outcomes', [], 'frequencies', []);
            continue;
        end
        
        % Extract the columns corresponding to the current subset
        selectedColumns = binaryVectors(:, subset);
        subsetSize = sum(subset);
        
        % Convert each row to a unique integer for easy counting
        powers = 2.^(subsetSize-1:-1:0);
        outcomeIndices = selectedColumns * powers';
        
        % Count the occurrences of each outcome
        [uniqueOutcomes, ~, outcomeIdx] = unique(outcomeIndices);
        frequencies = histcounts(outcomeIdx, 1:length(uniqueOutcomes)+1);
        
        % Convert unique outcomes back to binary using bitget
        binaryOutcomes = zeros(length(uniqueOutcomes), subsetSize);
        for j = 1:length(uniqueOutcomes)
            for k = 1:subsetSize
                binaryOutcomes(j, k) = bitget(uniqueOutcomes(j), subsetSize - k + 1);
            end
        end
        
        % Save results
        result{i+1} = struct('subset', find(subset), ...
                             'outcomes', binaryOutcomes, ...
                             'frequencies', frequencies, ...
                             'relativeFrequencies', frequencies/N);
    end
end

