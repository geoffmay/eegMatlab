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


robi003Icafilenames = [...
  {'/home/data/EEG/processed/Robi/robi003icaPost.mat'}, ...
  {'/home/data/EEG/processed/Robi/robi003icaPre.mat'}, ...
  ];
if(exist(robi003Icafilename, 'file'))
  load(robi003Icafilename);
else
  
  %filename = '/Users/Geoff/Documents/MATLAB/EEG/Coe Collection/Robi/ROBI_003/tx 1/630165006453007385.eegData';
  %        filename = '/home/data/EEG/data/ROBI/ROBI_003/outcome eyes open/630230539149591228.eegData';
  dataLength = NaN(length(filenames),1);
  allData = [];
  for fileCounter = 1:length(filenames)
    filename = filenames{fileCounter};
    outputStart = strfind(filename, 'ROBI_');
    if(length(outputStart) < 1)
      error('invalid filename');
    end
    
    channelCount = 34;
    file = dir(filename);
    fileLength = file.bytes / 8;
    fileId = fopen(filename);
    
    contents = fread(fileId, fileLength, 'double');
    
    fclose(fileId);
    data = reshape(contents, channelCount, fileLength / channelCount)';
    data = data(:, 1:32);
    dataLength(fileCounter) = size(data,1);
    allData = [allData; data];
  end
  clear contents;
  clear data;
  [ica.icaSig, ica.mixingMatrix, ica.separatingMatrix] = fastica(allData', 'maxNumIterations', 150);
  ica.filenames = filenames;
  ica.dataLengths = dataLength;
  labels = antChannelLocs;
  ica.channelLabels = labels(1:32);

  save(robi003Icafilename, 'ica');
  
end

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
  end
end

dataCounter = 1;
for fileCounter = 1:length(ica.dataLengths)  
  channelPair.coherence = [];
  channelPair.label = '';
  dataEnd = dataCounter + ica.dataLengths(fileCounter) - 1;
  icaData = ica.icaSig(dataCounter:dataEnd, :);
  pairCounter = 1;
  for i = 1:(maxChannel)
    for j = i+1:(maxChannel)
      fprintf('%d, %d\n',i,j);
      channelPair.coherence = coherence(icaData(:,i), icaData(:,j));
      %debug
      close all;
      figure;
      x = 1:size(channelPair.coherence,1);
      x = x / 128;
      plot(x, channelPair.coherence);
      xlabel('time (seconds)');
      ylabel('coherence');
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
  outputPath =robi003Icafilenames{fileCounter};
  save(outputPath, 'channelPairs', 'ica');
  dataCounter = dataCounter + ica.dataLengths(fileCounter);
end


