function [ pairs ] = getRobiCompleterFiles( onlyVerum )
%GETROBICOMPLETERFILES Summary of this function goes here
%   Detailed explanation goes here

if(~exist('onlyVerum', 'var'))
  onlyVerum = 1;
end

allFiles = getRobiDataFiles();
baseI = cellfun(@length, strfind(allFiles, 'baseline eyes open')) > 0;
exitI = cellfun(@length, strfind(allFiles, 'outcome eyes open')) > 0;
artifactI = cellfun(@length, strfind(allFiles, 'artifact')) > 0;

base = allFiles(baseI & ~artifactI);
exit = allFiles(exitI & ~artifactI);

pairs = cell(0);
for i = 1:length(exit)
  subject = exit{i};
  ind = strfind(subject,'ROBI_');
  subject = subject(ind:ind+7);
  baseCross = find(cellfun(@length, strfind(base, subject)) > 0);
  pairs(i, 1) = base(baseCross);
  pairs(i, 2) = exit(i);
end

if(onlyVerum)
  drop = [8,9];
  remove = logical(zeros(size(pairs,1),1));
  for i = 1:length(drop)
    filter = sprintf('%03d', drop(i));
    hit = find(cellfun(@length, strfind(pairs(:,1), filter)) > 0);
    remove(hit) = 1;
  end
  pairs(remove,:) = [];
end

end

