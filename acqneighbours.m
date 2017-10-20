function[PE_index]=acqneighbours(s,b,Npe)
%s is the sampling factor 
%b is the number  of acquired neighbours IN PE DIRECTION
%Npe is the number of phase encoding lines
PE_index=[];
index = 1;
while (1 - s + s * (max(1, 1 + index - floor(b / 2)) + b - 1) <= Npe)
    PE_index(index * s - s + 1:index * s, :) = repmat(1 - s + s * [max(1, 1 + index - floor(b / 2)) : max(1, 1 + index - floor(b / 2)) + b - 1] , [s  1]);
    index = index + 1;
end
for index = size(PE_index, 1) + 1:Npe
    PE_index(index, :) = PE_index(index - 1, :);
    PE_index(index, end) = Npe + 1;
end