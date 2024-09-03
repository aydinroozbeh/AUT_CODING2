% This code calculates the message from parity check node to check node
% M = current matrix of M's
% E = current matrix of E's
% i = related variable node
% j = related parity check node
    function [E_out] = cal_e(M_in)
    
    % Extracting dimensions
    n = size(M_in,2);
    k = n - size(M_in,1);
    E_out = zeros(size(M_in));
    
    row_mul = ones(n-k,1);
    M_tan = tanh(M_in/2);
    
    for j=1:1:(n-k)
        for i=1:1:n
            if(M_in(j,i)~=0)
                row_mul(j)=row_mul(j) * M_tan(j,i);
            end
        end
    end
    
    for j=1:1:(n-k)
        for i=1:1:n
            if(M_in(j,i)~=0)
                temp = row_mul(j)/ M_tan(j,i);
                E_out(j,i) = log((1+temp)/(1-temp));
            end
        end
    end