% Calculate the number of frame error rates

function num = cal_frame_error(A,B)

n = size(A,1);
num=0;

for i=1:1:n
    if( ~isequal(A(i,:) , B(i,:)) )
        num=num+1;
    end
end