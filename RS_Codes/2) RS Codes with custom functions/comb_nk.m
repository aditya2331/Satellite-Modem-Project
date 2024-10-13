% Function to generate error vectors
function result = comb_nk(n, k)
    result = [];
    
    % Total number of possible combinations
    num_combinations = numel(n)^k;
    
    % Generate all combinations
    for i = 1:num_combinations
        % Convert index to combination
        index = i - 1;
        combination = zeros(1, k);
        for j = k:-1:1
            combination(j) = n(mod(index, numel(n)) + 1);
            index = floor(index / numel(n));
        end
        result = [result; combination];
    end
end