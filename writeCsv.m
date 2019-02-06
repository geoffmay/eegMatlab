function writeCsv(filename, tab)

fileId = fopen(filename, 'w');

for i = 1:length(tab.Properties.VariableNames)
   content =  tab.Properties.VariableNames{i};
   fprintf(fileId, escape(content));
   if(i < length(tab.Properties.VariableNames))
       fprintf(fileId, ',');
   end
end

for i = 1:size(tab, 1)
    fprintf(fileId, sprintf('\n'));
    for j = 1:size(tab, 2)
        content = tab{i,j};
        if(iscell(content))
            content = content{1};
        end
        if(isnumeric(content))
            content = num2str(content);
            if(isnumeric(content))
                len = length(content);
            end
        end
        fprintf(fileId, escape(content));
        if(j < size(tab,2))
            fprintf(fileId, ',');
        end
    end
end

fclose(fileId);

end

function output = escape(input)

%debug
a = whos('input');
fprintf('\n %s: %s', a.class, input);
%end debug

output = input;
if(contains(input, '"') || contains(input, sprintf('\n')) || contains(input, ','))
    output = strrep(output, '"', '""');
    output = ['"' output '"'];
end
output = strrep(output, '%', '%%');

end