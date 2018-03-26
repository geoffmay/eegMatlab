function [ number ] = parseInt( text )
%PARSEINT Summary of this function goes here
%   Detailed explanation goes here

number = 0;
cursor = 1;
textLength = length(text);
badInput = false;
while(cursor <= textLength && ~badInput)
    c = text(cursor);
    if('0' <= c && c <= '9')
        number = number * 10 + (c - '0');
    else
        badInput = true;
    end
    cursor = cursor + 1;
end
if(badInput && cursor == 2)
    error('invalid input text');
end

end

