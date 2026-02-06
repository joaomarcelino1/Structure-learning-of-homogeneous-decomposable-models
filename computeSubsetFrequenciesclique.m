function result = computeSubsetFrequenciesclique(binaryVectors, kappa)
    % Input: Dataset D
    % Output: result - A cell array where each cell corresponds to a subset of coordinates
    %         and contains frequency information for that subset
    [N, n] = size(binaryVectors);
    
    % Validate kappa
    if kappa < 0 || kappa > n
        error('kappa must be between 0 and n');
    end
    
    % Count total number of valid subsets (with size <= kappa)
    totalSubsets = 0;
    for k = 0:kappa
        totalSubsets = totalSubsets + nchoosek(n, k);
    end
    
    result = cell(totalSubsets, 1);
    resultIndex = 1;
    
    for subsetSize = 0:kappa
        % Generate all combinations of the current size
        if subsetSize == 0
            % Empty set
            result{resultIndex} = struct('subset', [], 'outcomes', [], 'frequencies', []);
            resultIndex = resultIndex + 1;
        else
            % Generate all combinations of subsetSize elements from 1:n
            combinations = nchoosek(1:n, subsetSize);
            
            for combIdx = 1:size(combinations, 1)
                subset = combinations(combIdx, :);
                
                selectedColumns = binaryVectors(:, subset);
                
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
                result{resultIndex} = struct('subset', subset, ...
                                           'outcomes', binaryOutcomes, ...
                                           'frequencies', frequencies, ...
                                           'relativeFrequencies', frequencies/N);
                resultIndex = resultIndex + 1;
            end
        end
    end
end