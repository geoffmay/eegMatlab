
clear;
text = 'what up what up what up home slice';
nGramSize = 2;


a = containers.Map();


spaces = strfind(text, ' ');
spaces = [0 spaces length(text)+1];

spaceIndex = 1;
endSpace = spaceIndex + nGramSize;
while(endSpace <= length(spaces))
    if(mod(endSpace,1000)==0)
        fprintf('\n%f', endSpace / length(spaces));
    end
    nGram = text(spaces(spaceIndex)+1:spaces(endSpace)-1);
    if(isKey(a, nGram))
        value = values(a, {nGram});
        value = value{1};
        a = [a; containers.Map({nGram}, value + 1)];
    else
        if(a.Count == 0)
            a = containers.Map({nGram}, 1);
        else
            a = [a; containers.Map({nGram}, 1)];
        end
    end
    spaceIndex = spaceIndex + 1;
    endSpace = spaceIndex + nGramSize;
end