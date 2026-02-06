function sortedResult = sortSubsetResultsclique(result,nvar)
% Sort the result: [], [1], [2], ..., [1,2], [1,3], ...
% Input: result - output of computeSubsetFrequencies
% Output: sortedResult 

numCells = length(result);

n = nvar;

keyMatrix = zeros(numCells, n+1);


originalIndices = 1:numCells;

for i = 1:numCells
    
    subset = result{i}.subset;

    keyMatrix(i, 1) = length(subset);
    
    if ~isempty(subset)
        maxIdx = min(length(subset), n);
        keyMatrix(i, 2:(maxIdx+1)) = subset(1:maxIdx);
    end
end

[~, sortOrder] = sortrows(keyMatrix);

sortedResult = cell(numCells, 1);
for i = 1:numCells
    sortedResult{i} = result{originalIndices(sortOrder(i))};
end
end