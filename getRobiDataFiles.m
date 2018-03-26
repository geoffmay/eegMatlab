
function [paths, files, folders] = getRobiDataFiles(filter)
debug = true;
files = cell(0);
paths = cell(0);
folders = cell(0);
folder = '/home/data/EEG/data/ROBI';
subdirs = dir(folder);
for i = 1:length(subdirs)
  if(~strcmp(subdirs(i).name, '.') & ~strcmp(subdirs(i).name, '..'))
    subdir = fullfile(folder, subdirs(i).name);
    proceed = true;
    if(exist('filter','var'))
      proceed = false;
      if(length(strfind(subdir,filter)) > 0)
        proceed = true;
      end
    end
    if(proceed)
      subdirs2 = dir(subdir);
      hitCount = 0;
      for j = 1:length(subdirs2)
        if(~strcmp(subdirs2(j).name, '.') & ~strcmp(subdirs2(j).name, '..'))
          subdir2 = fullfile(subdir, subdirs2(j).name);
          if(length(strfind(subdir2, 'tx')) > 0 | length(strfind(subdir2, 'eyes')) > 0)
            if(length(strfind(subdir2, 'impedance')) == 0)
              theseFiles = dir(subdir2);
              for k = 1:length(theseFiles)
                if(strfind(theseFiles(k).name, '.eegData'))
                  if(theseFiles(k).name(1) ~= '.')
                    files{end+1} = theseFiles(k).name;
                    folders{end+1} = subdir2;
                    paths{end+1} = fullfile(subdir2, theseFiles(k).name);
                    hitCount = hitCount + 1;
                  end
                end
              end
            end
          end
        end
      end
      fprintf('\n%s: %d',subdir, hitCount);
    end
  end
end
fprintf('\ntotal: %d\n',length(files));

end