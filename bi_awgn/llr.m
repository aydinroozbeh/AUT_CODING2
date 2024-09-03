% Calculate the log likelihood ratio for two possilbe channels
% LLR is calcualted bit by bit - BPSK is assumed

function [out] = llr(y, SNR_db)

% calculating the absolute value
SNR_temp = 10^(SNR_db/10);

out = 4 * SNR_temp  * y;


end


