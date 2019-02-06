function [ output_args ] = saveRobiDataFile( filename, data )
%SAVEROBIDATAFILE Summary of this function goes here
%   Detailed explanation goes here

fileId = fopen(filename, 'wb');

if(size(data,1) < size(data,2))
  data = data';
end

for i = 1:size(data,1)
  for j = 1:size(data,2)
    fwrite(fileId, data(i,j), 'double');
  end
end

fclose(fileId);

end

