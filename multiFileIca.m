
maxRam = 10 * 1024 * 1024 * 1024;

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

toProcess = [pairs(:,1); pairs(:,2)]; 
masterIcaCoh = deriveIcaCoherenceMatrix(toProcess, 32);
outPath = '/home/data/EEG/processed/masterIcaCoh.mat';
save(outPath, 'masterIcaCoh', '-v7.3');
% 
% totalSize = 0;
% for i = 1:size(pairs,1)
%   for j = 1:size(pairs,2)
%     fileInfo = dir(pairs{i,j});
%     totalSize = totalSize + fileInfo.bytes; 
%   end
% end
% 
% if(totalSize > maxRam)
%   error('files exceed ram limits; size reduction not yet implemented')
% end
% 
% fprintf('\n%s: reading data', char(datetime));
% dataCounter = 1;
% allData = [];
% for i = 1:size(pairs,1)
%   for j = 1:size(pairs,2)
%     filename = pairs{i,j};
%     data = loadRobiDataFile(filename);
%     dataLength = size(data,1);
%     allData(dataCounter:dataCounter+dataLength-1,:) = data;
%     dataCounter = dataCounter + dataLength;
%   end
% end
% 
% fprintf('\n%s: performing ica', char(datetime));
% %todo:
% icaCoh = deriveIcaCoherenceMatrix();
