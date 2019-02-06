%  eegDataFilename = '/home/data/EEG/data/RobiPilot/RA/baseline eyes open/63055198303000001.eegData';


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


for fileCounter = 1:size(pairs, 1)
  
  
  
  for permCounter = 1:3
    
    if(permCounter == 1)
      eegDataFilename = pairs(fileCounter,:);
      desc = 'combined';
    elseif(permCounter == 2)
      eegDataFilename = pairs(fileCounter,1);
      desc = 'pre'
    elseif(permCounter == 3)
      eegDataFilename = pairs(fileCounter,2);
      desc = 'post'
    end
    
    slashes = strfind(eegDataFilename{1}, '/');
    subjectName = [eegDataFilename{1}(slashes(end-2)+1:slashes(end-1)-1) desc];
    outputFilename = ['/home/data/EEG/processed/Robi/icaCohPca/' subjectName '_ica.csv'];
    pcaFilename = ['/home/data/EEG/processed/Robi/icaCohPca/' subjectName  '_pca.csv'];
    
    
    % saveRobiDataFile(piece1, data1);
    
    fprintf('\n(%s) starting processing', char(datetime));
    tic;
    processIcaCohPca1(eegDataFilename);
    toc;
    
    
  end
end

