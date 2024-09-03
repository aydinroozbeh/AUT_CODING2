% convert LLR to bits

function [code] = get_code(LLR)

n=size(LLR,2);
code = zeros(1,n);

for i=1:1:n
    if(LLR(i)>0)
        code(i)=0;
    else
        code(i)=1;
    end
end