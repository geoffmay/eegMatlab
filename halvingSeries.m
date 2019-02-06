function [output] = halvingSeries(series)
%HALVINGSERIES Summary of this function goes here
%   Detailed explanation goes here

output = NaN(size(series));
taken = false(size(series));
output(1) = series(1);
taken(1) = true;
output(2) = series(end);
taken(end) = true;
writeIndex = 3;
divider = 2;
while(divider < length(series))
    increment = length(series) / divider;
    for i = increment:(increment):(length(series))
        ind = floor(i);
        if(~taken(ind))
            output(writeIndex) = series(ind);
            taken(ind) = true;
            writeIndex = writeIndex + 1;
        end
    end
    divider = divider * 2;
end
for i = 1:length(taken)
    if(~taken(i))
        output(writeIndex) = series(i);
        writeIndex = writeIndex + 1;
    end
end

end

