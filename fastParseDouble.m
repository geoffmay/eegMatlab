function number = fastParseDouble (text)

coefficient = 1;
number = 0;
decimalPlace = 1;

for i = 1:length(text)
    c = text(i);
    
    if('0' <= c && c <= '9')
        if(decimalPlace == 1)
        number = number * 10 + (c - '0');
        else
            number = number + decimalPlace * (c - '0');
            decimalPlace = decimalPlace * 0.1;
        end
    elseif(c == '-')
        coefficient = -1;
    elseif(c == '.')
        decimalPlace = 0.1;
    end
end
number = number * coefficient;



