
function sortedResult = sortSubsetResults(result)
% Sort as follows: [], [1], [2], ..., [1,2], [1,3], ...
% Input: result - output from computeSubsetFrequencies
% Output: sortedResult

numCells = length(result);

n = log2(numCells);

keyMatrix = zeros(numCells, n+1);

originalIndices = 1:numCells;

for i = 1:numCells
    subset = result{i}.subset;
    
    keyMatrix(i, 1) = length(subset);
    % Empty sets, columns are zero 
    if ~isempty(subset)
        % Do not exceed matrix
        maxIdx = min(length(subset), n);
        keyMatrix(i, 2:(maxIdx+1)) = subset(1:maxIdx);
    end
end

% Sort
[~, sortOrder] = sortrows(keyMatrix);
sortedResult = cell(numCells, 1);
for i = 1:numCells
    sortedResult{i} = result{originalIndices(sortOrder(i))};
end
end
