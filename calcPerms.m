% Calculate all permutations of a set
function Perms = calcPerms(setS)
    %INPUT: setS
    %OUTPUT: Vector with all permutations of setS
    maxSize = length(setS) * 24;
    
    % Pre-locate 
    Perms = cell(1, maxSize);
    totalIdx = 0;
    

    for i = 1:length(setS)

        currentVec = setS{i};
        
        % Calculate all perms
        allPerms = perms(currentVec);
        
        % CConvert to string format
        permStrings = cellfun(@(x) num2str(x(:)'), mat2cell(allPerms, ones(size(allPerms,1),1), length(currentVec)), 'UniformOutput', false);
        
        % Remove duplicate ones
        [uniquePermStrings, uniqueIdx] = unique(permStrings);
        uniquePerms = allPerms(uniqueIdx, :);
        
        % Add unique permutations
        numUnique = length(uniquePermStrings);
        Perms(totalIdx+1:totalIdx+numUnique) = mat2cell(uniquePerms, ones(numUnique,1), length(currentVec));
        totalIdx = totalIdx + numUnique;
    end
    
    % Remove unused pre-located memory
    Perms = Perms(1:totalIdx);
end