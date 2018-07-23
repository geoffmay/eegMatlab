function output = makeTableNumeric(input)

debug = false;
if(debug)
  debugCounter = 1;
end

for columnCounter = 1:size(input,2)
  column = input{:, columnCounter};
  if(~isnumeric(column))
    converted = NaN(size(column));
    goodParse = true;
    for i = 1:length(column)
      item = column{i};
      if(strcmp(item, 'null'))
      elseif(strcmp(item, ''))
      elseif(strcmp(item, 'NaN'))
      else
        number = str2double(item);
        if(~isnan(number) & (length(number) > 0))
          converted(i) = number;
        else
          goodParse = false;
        end
      end
    end
    if(goodParse)
      columnHeader = input.Properties.VariableNames{columnCounter};
      if(columnCounter > 1)
        preTable = input(:, 1:columnCounter-1);
      else
        preTable = table;
      end
      if(columnCounter < size(input, 2))
        postTable = input(:, columnCounter+1:end);
      else
        postTable = table;
      end
      inTable = table(converted);
      inTable.Properties.VariableNames{1} = columnHeader;
      input = [preTable, inTable, postTable];
      if(debug)
        fprintf('\n(%d) replaced column %d', columnCounter, debugCounter);
        pre{debugCounter} = column;
        post{debugCounter} = inTable;
        debugCounter = debugCounter + 1;
      end
    end
  end
end

output = input;

%column = input{:, 'CRIS_Baseline_csvCRIS_Extent_1_4'};







