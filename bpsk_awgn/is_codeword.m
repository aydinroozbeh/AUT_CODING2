function [stat] = is_codeword(pcm , y)

n=size(pcm,2);
r=size(pcm,1);

if(isequal( mod(pcm*y',2) , zeros(r,1)))
    stat=1;
else
    stat=0;
end