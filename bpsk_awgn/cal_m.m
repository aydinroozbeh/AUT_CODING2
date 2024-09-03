% Calcualtes the message from variable node to parity check node
% M = current matrix of M's
% E = current matrix of E's
% i = related variable node
% j = related parity check node

function [M_out] = cal_m(E_in , r)

E_sum = sum(E_in);
M_out = zeros(size(E_in));
 
% Extracting dimensions
n = size(E_in,2);
k = n - size(E_in,1);

for j=1:1:(n-k)
    for i=1:1:n
        if(E_in(j,i)~=0)
            M_out(j,i)=E_sum(i)-E_in(j,i)+r(i);
        end
    end
end