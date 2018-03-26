function [ output ] = transposeTable( input )
%TRANSPOSE Summary of this function goes here
%   Detailed explanation goes here

output = cell(size(input,2),size(input,1));
for(i = 1:size(input,1))
    for(j = 1:size(input,2))
        output{j,i}=input(i,j);
    end
end

end

