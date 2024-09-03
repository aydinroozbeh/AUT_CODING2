% Sum-product decoding algorithm
% Variable definition:
% Decoded: the decoded codeword after running the algorithm
% Status: states whether the decoding has been successful
% max_iter : Maximum number of iterations
% r : received bits (the noisy input)
% the SNR of the channel(db), it is used to find the LLR values
% -----
function [decoded] = decoder_awgn(max_iter , pcm , y , SNR)

% Extracting the rate
n = size(pcm , 2);
k = n - size(pcm,1);

% Number of iterations done up until now:
iter = 0;

% initializing LogLikelihood ratios based on received signal values
r=zeros(1,n);
for i=1:1:n
    r(i) = llr(y(i), SNR);
end

% Message Matrices - Have the same sparsity pattern as that of the PCM
M = zeros(size(pcm));
E = M;

% Current LLR

% Current Estimation of the received code word
temp = y;

% Initializing Matrix M based on received signal strength:
for j=1:1:(n-k)
    for i=1:1:n
        if(pcm(j,i)==1)
            M(j,i)=r(i);
        end
    end
end
    
while( (iter <= max_iter) && ~is_codeword(pcm,temp) )
    if(iter~=0)
        M = cal_m(E,r);     % Updating the 'M' Matrix
    end
        E = cal_e(M);       % Updating the 'E' Matrix
        L = sum(E) + r;     % Calculating new LLRs
        temp = get_code(L); % Getting the bits corresponding to LLRs
        
    iter = iter +1;


end

decoded = temp;
