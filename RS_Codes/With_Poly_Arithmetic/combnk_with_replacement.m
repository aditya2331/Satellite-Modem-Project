% Function to generate error vectors
function result = combnk_with_replacement(v, k)
    result = [];
    
    % Total number of possible combinations
    num_combinations = numel(v)^k;
    
    % Generate all combinations
    for i = 1:num_combinations
        % Convert index to combination
        index = i - 1;
        combination = zeros(1, k);
        for j = k:-1:1
            combination(j) = v(mod(index, numel(v)) + 1);
            index = floor(index / numel(v));
        end
        result = [result; combination];
    end
end