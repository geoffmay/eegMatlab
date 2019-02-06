function [tab] = insertBlankRow(tab, rowIndex)
%INSERTBLANKROW Summary of this function goes here
%   Detailed explanation goes here

cursor = size(tab, 1);
while(cursor > rowIndex)
    tab(cursor + 1, :) = tab(cursor, :);
    cursor = cursor - 1;
end

for i = 1:size(tab, 2)
    thisCell = tab{cursor, i};
    if(isnumeric(thisCell))
        thisCell = NaN;
    elseif(iscell(thisCell))
        thisCell = cell(1);
    end
    tab{cursor,i} = thisCell;
end

end

