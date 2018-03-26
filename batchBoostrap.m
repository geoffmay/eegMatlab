createDummy = 1;

%get dyads of pre- and post-treatment files
files = getRobiDataFiles();
postIndexes = find(cellfun(@length, strfind(files, 'outcome')));
preIndexes = find(cellfun(@length, strfind(files, 'baseline')));
outputFilenames = cell(0);
openFilenames = cell(0);
closedFilenames = cell(0);
for i = 1:length(postIndexes)
  file = files{postIndexes(i)};
  target = strrep(file, 'outcome', 'baseline');
  target = target(1:strfind(target,' eyes ') + 8);
  targetIndexes = find(cellfun(@length, strfind(files, target)));
  info = file(strfind(file,'ROBI_'):end);
  slashes = strfind(info,'/');
  info(slashes(end):end) = [];
  info = strrep(info, '/outcome ','_');
  info = sprintf('/home/data/EEG/processed/Robi/coherenceBootstrapNoAutocorrelation/%s bootstrap.mat',info);
  if(length(strfind(info, 'closed'))==0)    
    openFilenames(size(openFilenames,1)+1, 1) = files(targetIndexes);
    openFilenames(size(openFilenames,1), 2) = files(postIndexes(i));
    openFilenames{size(openFilenames,1), 3} = info;
  else
    closedFilenames(size(closedFilenames,1)+1, 1) = files(targetIndexes);
    closedFilenames(size(closedFilenames,1), 2) = files(postIndexes(i));
    closedFilenames{size(closedFilenames,1), 3} = info;
  end
end
parameters = [openFilenames;closedFilenames];

for i = 1:length(parameters)
  file1 = parameters{i,1};
  file2 = parameters{i,2};
  output = parameters{i,3};
  fprintf('\n%s', output);
  if(~exist(output,'file'))
  if(createDummy)
    timeStarted = char(datetime);
    save(output, 'timeStarted');
  end
  data = bootstrap(file1, file2);
  save(output, 'data');
  end
  fprintf('\n');
end