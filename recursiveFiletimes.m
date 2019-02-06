function [ fileTimes ] = recursiveFiletimes( folder, fileTimes )
%RECURSIVEFILETIMES Summary of this function goes here
%   Detailed explanation goes here

contents = dir(folder);
for i = 1:length(contents)
  fullPath = fullfile(folder, contents(i).name);
  if(contents(i).isdir)
    if(~strcmp(contents(i).name, '.git') && ...
        ~strcmp(contents(i).name, '.') && ...
        ~strcmp(contents(i).name, '..'))
      fprintf('\n%s', fullPath);
      if(exist('fileTimes','var'))
        fileTimes = recursiveFiletimes(fullPath, fileTimes);
      else
        fileTimes = recursiveFiletimes(fullPath);        
      end
    end
  else
    output = sprintf('%f %s', contents(i).datenum, fullfile(folder, contents(i).name));
    fprintf('\n%s', output);
    if(~exist('fileTimes', 'var'))
      fileTimes = {output};      
    else
      fileTimes{end+1}=output;
    end
  end
end




end

