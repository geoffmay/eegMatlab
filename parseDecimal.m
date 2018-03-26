function [ sum ] = parseDecimal( text )
%PARSEDECIMAL Summary of this function goes here
%   Detailed explanation goes here

aspirationalPrecision = 0;
debug = 0;
if(text(1) == '-')
    sign = -1;
    startIndex = 2;
else
    sign = 1;
    startIndex = 1;
end
sum = 0;
pastDecimal = false;
decimalPlace = 1;
i = startIndex;
powerOfTen = 0;
while(i <= length(text))
    if(text(i) >= '0' && text(i) <= '9')
        number = text(i) - '0';
        if(~pastDecimal)
            sum = sum * 10 + number;
        else
            sum = sum + number * decimalPlace;
            decimalPlace = decimalPlace * .1;
            powerOfTen = powerOfTen -1;
        end
    elseif(text(i) == '.')
        pastDecimal = true;
        decimalPlace = .1;
    elseif(text(i) == 'I')
        sum = inf;
        break;
    else
        sum = NaN;
        break;
    end
    i = i + 1;
end
if(pastDecimal && aspirationalPrecision)
    hex = num2hex(sum);
    hex1 = hex(1:end-1);
    %todo: need to search a larger space (not sure how much larger,
    %possibly sweep based on finding a match very early or very late
    for i = 0:15
        c = i + '0';
        if(i > 9)
            c = i + 'a' - 10;
        end
        hex2 = [hex1 c];
        num2 = hex2num(hex2);
        text2 = sprintf('%.100f', num2);
        match = strfind(text2, text);
        if(match)
            if(exist('lastText','var'))
                cutoff = length(text) - 1;
                trunc1 = lastText(cutoff:end);
                trunc2 = text2(cutoff:end);
                target = text(cutoff:end);
                precision = 40;
                numTarget = parseDecimal(target) * power(10, precision - 2);
                resid1 = parseDecimal(trunc1(1:precision));
                resid2 = parseDecimal(trunc2(1:precision));
                diff1 = abs(resid1 - numTarget);
                diff2 = abs(resid2 - numTarget);
                if(diff1 > diff2)
                    sum = num2;
                else
                    sum = lastNum;
                end
                break;
            else
                sum = num2;
            end
        end
        lastHex = hex2;
        lastNum = num2;
        lastText = text2;
        if(debug)
            fprintf('\n%d %.100f', length(match), hex2num(hex2));
        end
    end
    if(~match)
        sum = lastNum;
    end
end

sum = sum * sign;
end



