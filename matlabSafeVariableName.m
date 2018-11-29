function [ variableName ] = matlabSafeVariableName( variableName )
%MATLABSAFEVARIABLENAME replace unsafe characters
%   Replaces [#] [.] [$] [ ] respectively with 
%   [NUM] [POINT] [DOLLAR] [_]
%             variableName = strrep(variableName, '#', 'NUM');
%             variableName = strrep(variableName, '.', 'POINT');
%             variableName = strrep(variableName, '$', 'DOLLAR');
            variableName = strrep(variableName, '#', '_');
            variableName = strrep(variableName, '.', '_');
            variableName = strrep(variableName, '$', '_');
            variableName = strrep(variableName, ' ', '_');
            variableName = strrep(variableName, '-', '_');
            variableName = strrep(variableName, '>', '_');
            variableName = strrep(variableName, '<', '_');
            variableName = strrep(variableName, '/', '_');
            variableName = strrep(variableName, '\', '_');
            variableName = strrep(variableName, '(', '_');
            variableName = strrep(variableName, ')', '_');
            variableName = strrep(variableName, ',', '');
            
            if(length(variableName > 0))
            unsafestarts = {'_','0','1','2','3','4','5','6','7','8','9'};
            if(sum(strcmp(unsafestarts, variableName(1))))
                variableName = strcat('a', variableName);
            end
            else
                variableName = 'a';
            end
            if(length(variableName) > namelengthmax)
                variableName = variableName(1:namelengthmax);
            end
end