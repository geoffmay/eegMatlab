function mat = loadRawMatrix(filename)
%LOADRAWMATRIX Summary of this function goes here
%   Detailed explanation goes here


fileId = fopen(filename, 'r');
rowCount = fread(fileId, 1, 'double');
colCount = fread(fileId, 1, 'double');
mat = NaN(rowCount, colCount);

for i = 1:colCount
    for j = 1:rowCount
        mat(j,i) = fread(fileId, 1, 'double');
    end
end

fclose(fileId);


end

