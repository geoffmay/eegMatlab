function cleanTemp()

rootFolder = '/tmp';
recursiveDelete(rootFolder);


function recursiveDelete(folder)
fprintf('\n%s', folder);
contents = dir(folder);
for i = 3:length(contents)
    item = contents(i);
    path = fullfile(folder, item.name);
    if(item.isdir)
        recursiveDelete(path);
        try
          rmdir(path);
        catch ex
          %nothing for now
        end
    else
        delete(path);
    end
end

