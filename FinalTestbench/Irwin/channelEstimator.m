%% Channel Estimator Function

function h_est = channelEstimator(rx_sampled,pilot_matrix)
    % rx_sampled(column vec): inital 26 bits must be the convolution of channel with
    % the pilot sequence

    % pilot_matrix : inverse of toeplitz matrix creaed from the best
    % preamble sequence
    rx_sampled = transpose(rx_sampled);
    h_est = transpose(pilot_matrix*transpose(rx_sampled(1,13:26))); % 
    h_est = transpose(h_est);
end



%% Supporting Function For Channel Estimator

% function to find the pseudo
function [Xp,invertible] = pseudoInvOfMat(seq)
    % Input : seq (row vec)
    % Output: Xp  (matrix)
    % Output: invertible (bool)
    % seq   : preamble sequence of length 26, elements should be from [0 1]

    % converting [0 1] sequence to [-1 +1] sequence
    for i = 1:length(seq)
        if(seq(1,i)==0)
            seq(1,i) = -1;
        end
    end

    % construct the channel convolution matrix
    X = createToeplitzMatrix(seq);

    A = transpose(X)*X;
    if(isempty(find(eig(A)==0,1)))
        invertible = true;
        X_inv = pinv(X);
    else
        invertible = false;
        X_inv = zeros(13,14);
    end
    Xp = X_inv; 
end

function T = createToeplitzMatrix(seq)
    % creates the 14x13 matrix from the preamble sequence
    T = zeros([14,13]);
    for i = 1:14
        T(i,:) = flip(seq(1,i:12+i));
    end 
end



