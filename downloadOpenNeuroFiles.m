url = 'https://openneuro.org/datasets/ds001512/versions/1.0.1';
%destRoot = 'C:\Users\Neuro\Downloads\GhermanPhiliastides\';
destRoot = 'C:\Users\Neuro\Documents\MATLAB\data\GhermanPhilastides\1.0.1';
sourceRoot = 'https://openneuro.org/crn/datasets/ds001512/snapshots/1.0.0/files/';


listFile = 'C:\Users\Neuro\Documents\files.txt';
fileId = fopen(listFile);
text = fscanf(fileId, '%c');
fclose(fileId);

linkMarker = 'href="';
linkStarts = strfind(text, linkMarker);
quotes = strfind(text, '"');

for i = 1:length(linkStarts)
    linkStart = linkStarts(i) + length(linkMarker);
    linkEnd = quotes(min(find(quotes > linkStart))) - 1;
    links{i} = text(linkStart:linkEnd);
end

isDisplay = cellfun(@length, strfind(links, 'file-display')) > 0;
fileLinks = links(~isDisplay);



for i = 1:length(fileLinks)
    link = fileLinks{i};
    sourceShort = link(length(sourceRoot)+1:end);
    destShort = strrep(sourceShort, ':', '\');
    destPath = fullfile(destRoot, destShort);
    if(~exist(destPath, 'file'))
        fprintf('\n%s: file %d of %d, (%s)', char(datetime), i, length(fileLinks), link);
        destParts = strsplit(destPath, '\');
        buildPath = destParts{1};
        for j = 2:length(destParts)-1
            buildPath = fullfile(buildPath, destParts{j});
            if(~exist(buildPath, 'file'))
                mkdir(buildPath);
            end
        end
        try
            websave(destPath, link);
        catch ex
            fprintf('\ncaught error: %s', ex.message);
            for stackCounter = 1:length(ex.stack)
                fprintf('\n  in: %s line %d', ex.stack(stackCounter).name, ex.stack(stackCounter).line);
            end
            delete(destPath);
        end
    end
end