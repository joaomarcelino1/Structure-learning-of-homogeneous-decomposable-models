% Generate all possible subsets of a set N
function subconjuntos = gerarsetS(N)
    %INPUT: set N
    %OUTPUT: all subsets of N
    total_subconjuntos = 2^N;
    subconjuntos = cell(1, total_subconjuntos);
    
    % For every binnary number from 0 to 2^N-1
    for i = 0:total_subconjuntos-1
        % Check elements
        bin = dec2bin(i, N);
        subset = find(bin == '1');
        
        % Get the numbers
        subset = N - subset + 1;
        subset = sort(subset); 
        
        subconjuntos{i+1} = subset;
    end
    
    % Sort by size, and then lexicographically.
    tamanhos = cellfun(@length, subconjuntos);
    [~, idx_ordenado] = sortrows([tamanhos', cell2mat(cellfun(@(x) [x, zeros(1, N-length(x))], subconjuntos, 'UniformOutput', false)')]);
    subconjuntos = subconjuntos(idx_ordenado);
end
