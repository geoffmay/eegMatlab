function saveRawMatrix(mat, filename)
%SAVERAWMATRIX Summary of this function goes here
%   Detailed explanation goes here
if(length(size(mat)) > 2)
    error('matrix has too many dimensions, must not have more than 2');
end

fileId = fopen(filename, 'w');
for i = 1:length(size(mat))
    fwrite(fileId, size(mat, i), 'double');
end

for i = 1:size(mat, 2)
    fwrite(fileId, mat(:,i), 'double');
end


fclose(fileId);

end

