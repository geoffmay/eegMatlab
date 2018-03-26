clear;

simOnly = false;
sampleRate = 2048;
epochLength = 1 * sampleRate;
i = 100;
startIndex = (i-1) * epochLength + 1;
endIndex = startIndex + epochLength - 1;

hiPassHz = 1;
loPassHz = 40;

filenames = [...
  {'/home/data/EEG/data/ROBI/ROBI_003/outcome eyes open/630230539149591228.eegData'}...
  {'/home/data/EEG/data/ROBI/ROBI_003/baseline eyes open/630158995692243270.eegData'}...
  ];

%filename = '/Users/Geoff/Documents/MATLAB/EEG/Coe Collection/Robi/ROBI_003/tx 1/630165006453007385.eegData';
%        filename = '/home/data/EEG/data/ROBI/ROBI_003/outcome eyes open/630230539149591228.eegData';


for fileCounter = 1:length(filenames)
  filename = filenames{fileCounter};
  outputFolder = '/home/data/EEG/processed/Robi/coherence/';
  outputStart = strfind(filename, 'ROBI_');
  if(length(outputStart) < 1)
    error('invalid filename');
  end
  outFile = filename(outputStart(1):end);
  outFile = sprintf('icaCoherence_%s', outFile);
  outFile = strrep(outFile,'/','_');
  outFile = strrep(outFile,'.eegData','.mat');
  outputPath = fullfile(outputFolder, outFile);

  channelCount = 34;
  file = dir(filename);
  fileLength = file.bytes / 8;
  fileId = fopen(filename);
  
  contents = fread(fileId, fileLength, 'double');
  
  fclose(fileId);
  data = reshape(contents, channelCount, fileLength / channelCount)';
  data = data(:, 1:32);
  clear contents;
  labels = antChannelLocs;
  
  [ica.icaSig, ica.mixingMatrix, ica.separatingMatrix] = fastica(data', 'maxNumIterations', 200);
  %  maxChannel = channelCount - 2;
  ica.icaSig = ica.icaSig';
  maxChannel = size(ica.icaSig, 2);  
  maxPairs = 0;
  indexTable = [];
  pairLabels = cell(0);
  for i = 1:(maxChannel)
    for j = i+1:(maxChannel)
      maxPairs = maxPairs + 1;
      indexTable(maxPairs, 1) = maxPairs;
      indexTable(maxPairs, 2) = i;
      indexTable(maxPairs, 3) = j;
%      pairLabels{maxPairs}= sprintf('%s-%s',labels{i},labels{j});
    end
  end
  
  channelPair.coherence = [];
  channelPair.label = '';  
  pairCounter = 1;
  for i = 1:(maxChannel)
    for j = i+1:(maxChannel)
      fprintf('%d, %d\n',i,j);
      channelPair.coherence = coherence(ica.icaSig(:,i), ica.icaSig(:,j));
      %debug
      close all;
      figure;
      plot(channelPair.coherence);
      legend('delta', 'theta', 'alpha', 'beta', 'hibeta');
      zoom xon;
      pan xon;
      %end debug
      if(any(any(isnan(channelPair.coherence))));
        fprintf('i = %d, j = %d\n', i, j);
        error('unexpected coherenceValue');
      end
      channelPair.label =  sprintf('ica%d-ica%d',i,j);
      fprintf('%d/%d\n',pairCounter,maxPairs);
      channelPairs(pairCounter) = channelPair;      
      pairCounter = pairCounter + 1;
    end
  end
  
  save(outputPath, 'channelPairs', 'filename', 'ica');
  
end
