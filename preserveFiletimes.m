root = fileparts(which('preserveFiletimes'));

files = recursiveFiletimes(root);

outputFile = fullfile(root, 'fileTimes.m');

fileId = fopen(outputFile, 'w');

for i = 1:length(files)
  line = files{i};
  line1 = strrep(line, root, '');
  fprintf(fileId, '\n%s', line1);
end

fclose(fileId);


