function output = eegCppAnalytics(EEG, dropChans)

refreshRate = 128;
fftWindowDuration = 1;
coherenceQueueDuration = 0.5;
freqs = '1 4 8 12 25 30';
%freqs = '1 5 9 13 26 30';
deleteTempFiles = 1;
debugPlot = 0;

tic;
oldFolder = pwd;

if(~exist('dropChans', 'var'))
  if(length(find(cellfun(@length, strfind({EEG.chanlocs.labels}, 'counter'))))>0)
    dropChans = {'M1', 'M2', 'CPz', 'counter'};
  else
    dropChans = {'EXG1', 'EXG2', 'EXG3', 'EXG4', 'EXG5', 'EXG6', 'EXG7', 'EXG8', 'Ana1', 'Ana3', 'Status', 'Cz'};
  end
end

if(exist('data', 'var'))
  EEG.data = data';
end

if(~exist('EEG', 'var'))
  error('must pass an EEG variable');
  data = loadRobiDataFile('/home/data/EEG/data/ROBI/ROBI_003/baseline eyes open/630158995692243270.eegData');
  EEG.data = data';
  EEG.srate = 2048;
  [~, EEG.chanlocs] = antChannelLocs;
  allChans = 'Fp1 Fpz Fp2 F7 F3 Fz F4 F8 FC5 FC1 FC2 FC6 M1 T7 C3 Cz C4 T8 M2 CP5 CP1 CP2 CP6 P7 P3 Pz P4 P8 POz O1 Oz O2 CPz counter';
  keepChans = 'Fp1 Fpz Fp2 F7 F3 Fz F4 F8 FC5 FC1 FC2 FC6 T7 C3 Cz C4 T8 CP5 CP1 CP2 CP6 P7 P3 Pz P4 P8 POz O1 Oz O2';
else
  chans = {EEG.chanlocs.labels};
  allChans = '';
  for i = 1:length(chans)
    allChans = [allChans chans{i}];
    if(i < length(chans))
      allChans = [allChans ' '];
    end
  end
  
  keep = true(1, length(EEG.chanlocs));
  for i = 1:length(keep)
    for j = 1:length(dropChans)
      if(strcmp(EEG.chanlocs(i).labels, dropChans{j}))
        keep(i) = 0;
      end
    end
  end
  EEG.data(~keep, :) = [];
  
  chans(~keep) = [];
  keepChans = '';
  for i = 1:length(chans)
    keepChans = [keepChans chans{i}];
    if(i < length(chans))
      keepChans = [keepChans ' '];
    end
  end
  
end

if(exist('data', 'var'))
  EEG.data = data';
end


currentFolder = fileparts(which('eegCppAnalytics'));
executeFolder = fullfile(currentFolder, 'cppFft');
cd(executeFolder);
binaryPath = fullfile(currentFolder, 'cppFft', 'computePowerAndCoherence.bin');
%binaryPath2 = fullfile(currentFolder, 'cppFft', 'a');
configFile = fullfile(currentFolder, 'cppFft', 'config.txt');

mode = 'offline';
numberOfFiles = 1;
collisionCounter = 0;
inputFile = fullfile(currentFolder, 'cppFft', 'temp', sprintf('temp%d.eegData', collisionCounter));
while(exist(inputFile,'file'))
  collisionCounter = collisionCounter + 1;
  inputFile = fullfile(currentFolder, 'cppFft', 'temp', sprintf('temp%d.eegData', collisionCounter));
end
outputFile = fullfile(currentFolder, 'cppFft', 'temp', sprintf('temp%d.analyte', collisionCounter));
sampleRate = EEG.srate;

%write input config file
fileId = fopen(configFile, 'w');
fprintf(fileId, '%s\n', mode);
fprintf(fileId, '%f\n', numberOfFiles);
fprintf(fileId, '%s\n', inputFile);
fprintf(fileId, '%s\n', outputFile);
fprintf(fileId, '%f\n', sampleRate);
fprintf(fileId, '%f\n', refreshRate);
fprintf(fileId, '%f\n', fftWindowDuration);
fprintf(fileId, '%f\n', coherenceQueueDuration);
fprintf(fileId, '%s\n', allChans);
fprintf(fileId, '%s\n', keepChans);
fprintf(fileId, '%s\n', freqs);
fclose(fileId);

%write input data
fileIdBin = fopen(inputFile, 'w');
a = reshape(EEG.data, [1, numel(EEG.data)]);
fwrite(fileIdBin, a, 'double');
fclose(fileIdBin);
system(binaryPath);

%read output header
fileId = fopen([outputFile '.txt']);
header = fscanf(fileId, '%c');
lines = strfind(header, sprintf('\n'));
output.sRate = sscanf(header(1:(lines(1)-1)), 'sample rate: %f');
output.cohRate = sscanf(header((lines(1)+1):(lines(2)-1)), 'coherence refresh rate: %f');
interp = header((lines(2)+1):(lines(3)-1));
output.interpolation = interp(length('interpolation: ')+1:end);
output.measureCount = sscanf(header((lines(3)+1):(lines(4)-1)), 'measure count: %f');
measures = header((lines(4)+1):(lines(5)-1));
measures = measures(length('measurements (comma separated): ')+1:end);
output.measures = strsplit(measures, ',');
fclose(fileId);


output.data = loadBinaryMatrix(outputFile, output.measureCount);
wastedColumns = find(all(isnan(output.data),1));
incompleteRows = find(any(isnan(output.data),2));
output.data(incompleteRows,:) = [];

if(deleteTempFiles)
  delete([outputFile '.txt']);
  delete(outputFile);
  delete(inputFile);
end

cd(oldFolder);

if(debugPlot)
  folder = '/home/data/EEG/processed/Robi/coherence3';
  save(fullfile(folder, 'binaryVersion.mat'), 'output')
  
  comp.mat = load(fullfile(folder, 'matlabVersion.mat'));
  comp.bin = load(fullfile(folder, 'binaryVersion.mat'));
  
  mat.matrix = comp.mat.coh.matrix;
  mat.labels = comp.mat.coh.labels;
  bin.matrix = comp.bin.output.data;
  bin.labels = comp.bin.output.measures;
  
  %bin is the binary produced data, mat is matlab produced.
  matI = 2329;
  binI = 4;
  matI = 4;
  binI = 154;
  labelM = mat.labels{matI};
  labelB = bin.labels{binI};
  dataM = mat.matrix(:,matI);
  dataB = bin.matrix(:,binI);
  len = min(length(dataM), length(dataB));
  dataM = dataM(1:len);
  dataB = dataB(1:len);
  dataB1 = log10(dataB);
  x = bin.matrix(1:len,end);
  plot(x, [dataM, dataB]);
  legend({'mat', 'bin'});
  title(sprintf('%s vs. %s', labelM, labelB));
  veryHigh = dataB > 0.99;
  top = dataB == 1;
  pan xon;
  zoom xon;
  
  figure;
  [pcaB.COEFF, pcaB.SCORE, pcaB.LATENT, pcaB.TSQUARED, pcaB.EXPLAINED, pcaB.MU] = pca(bin.matrix(:,1:end-1));
  [pcaM.COEFF, pcaM.SCORE, pcaM.LATENT, pcaM.TSQUARED, pcaM.EXPLAINED, pcaM.MU] = pca(mat.matrix(:,1:end));
  fftB = fft(pcaB.COEFF(:,1));
  fftM = fft(pcaM.COEFF(:,1));
  fftx = (1:len) .* (EEG.srate/len);
  ind = 2:1000;
  plot(fftx(ind), abs(fftB(ind)));
  plot(fftx(ind), abs(fftM(ind)));
  
end

end




