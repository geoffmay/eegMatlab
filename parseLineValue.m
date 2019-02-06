function [value] = parseLineValue(text, key)
%PARSELINEVALUE Summary of this function goes here
%   Detailed explanation goes here
headerLines = strsplit(text, newline);
keyLine=headerLines{find(cellfun(@length, strfind(headerLines, key))>0)};
valueText = keyLine(strfind(keyLine, key) + length(key):end);
valueText = strrep(valueText, sprintf('\r'), '');
value = str2double(valueText);

if(isnan(value))
    value = valueText;
end

end

