%Tests granger causality and nipal PCA on EEG data.

doGranger = true;
doPCA = true;
PcaComponentCount = 32;

%get list of files
folder = '/Users/Geoff/Documents/MATLAB/EEG/rebuilt feedback'; %mac
folder = '/Volumes/My Book/zScoresComplete'; %mac ext hdd
folder = '/home/data/EEG/processed/Robi/fullZ'; %server
folder = '/media/HERMANOSBAC/data/EEG1' %external drive
files = dir(folder);
files = files(~[files.isdir]);

%get variable info from one of the notes files.
searchString = '.zScorenotes.txt';
searchHits = strfind({files.name}, searchString);
searchHits = find(cellfun(@length, searchHits));
fileId = fopen(fullfile(folder,files(searchHits(1)).name));
text = textscan(fileId, '%s');
text1 = text{1};
text1a = text1(42:49);
text1b = text1(50:end);
text2 = [text1b;text1a];
fclose(fileId);
for i = 1:length(text2)
  text = text2{i};
  text = text(1:end-1);
  text2{i} = text;
end
columnLabels = text2;

%remove the non-data files
remove = zeros(length(files),1);
for i = 1:length(files)
  filename = files(i).name;
  if(filename(1) == '.')
    remove(i) = 1;
  end
  if(length(filename) > 11)
    if(filename(end-15:end) == '.zScorenotes.txt')
      remove(i) = 1;
    end
    if(filename(end-19:end) == '.zScorenotes.txt.zip')
      remove(i) = 1;
    end
  end
end
files(find(remove)) = [];

%process all files
for i = 1:length(files)
  %execute only if not already done
  outputFolder = fullfile(folder, 'pca');
  path = fullfile(outputFolder, sprintf('%snipals.mat',files(i).name));
  if(~exist(path))
    
    %shape data matrix based on column headers from info file
    fileId = fopen(fullfile(folder, files(i).name));
    contents = fread(fileId, files(i).bytes/4, 'single');
    fclose(fileId);
    columns = length(columnLabels);
    rows = length(contents) / columns;
    contents = reshape(contents, columns, rows)';
    maxJ = size(contents, 2);
    tickSize = floor(maxJ/100);
    if(doGranger) %granger causality
      tic;
      correlationRs = NaN(1,maxJ);
      correlationPs = NaN(1,maxJ);
      grangerF = NaN(maxJ, 1);
      grangerCrit = NaN(maxJ, 1);
      grangerCausal = NaN(maxJ, 1);
      for j = 1:maxJ
        fprintf('.');
        if(mod(j,tickSize) == 0)
          fprintf('.');
          if(mod(j,tickSize * 10))
            fprintf('%d', j / tickSize);
          end
        end
        max_lag = 2;
        alpha = 0.05;
        x = contents(:,j);
        y = contents(:,maxJ-3); %feedback
        [correlationRs(j), correlationPs(j)] = corr(x,y);
        [grangerF(j), grangerCrit(j)] = granger_cause(x,y,alpha,max_lag);
        grangerCausal(j) = grangerF(j) > grangerCrit(j);
      end
      elapsed = toc;
    end %granger causality
    
    if(false) %diy labels
      
      freqLabels = [{'Delta1to4Hz'},{'Theta4to8Hz'},{'Alpha8to12Hz'},{'Beta12to25Hz'},{'HiBeta25to30Hz'}];
      chanLabels = [{'Fp1'},{'Fp2'},{'F3'},{'F4'},{'C3'},{'C4'},{'P3'},{'P4'},{'O1'},{'O2'},{'F7'},{'F8'},{'T3'},{'T4'},{'T5'},{'T6'},{'Fz'},{'Cz'},{'Pz'}];
      chanPair = cell(0);
      for j = 1:length(chanLabels)
        for k = j+1:length(chanLabels)
          chanPair{end+1} = sprintf('%s-%s', chanLabels{j}, chanLabels{k});
        end
      end
      
      columnCounter = 1;
      for j = 1:length(freqLabels)
        for k = 1:length(chanLabels)
          columnLabels{columnCounter, 1} = sprintf('%s-%s', chanLabels{k}, freqLabels{j});
          columnCounter = columnCounter + 1;
        end
      end
      for j = 1:length(freqLabels)
        for k = 1:length(chanPair)
          columnLabels{columnCounter, 1} = sprintf('%s-%s', chanPair{k}, freqLabels{j});
          columnCounter = columnCounter + 1;
        end
      end
    end
    
    if(doGranger)
      causalRatio = grangerF ./ grangerCrit;
      standardDev = std(contents,1)';
      meanValue = mean(contents,1)';
      normalizedCausalRatio = causalRatio ./ standardDev;
      tab = table(columnLabels, grangerF, grangerCrit, grangerCausal, causalRatio, meanValue, standardDev, normalizedCausalRatio);
      save(path, 'dataPca', 'tab', 'max_lag', 'alpha', 'elapsed');
    else
      contents = contents(:, 1:end-10);
      tic;
      fprintf('\n\ndoing pca on file %d of %d', i, length(files));
      dataPca.maxIterations = power(10,6);
      dataPca.desiredComponents = 20;
      [dataPca.T, dataPca.P, dataPca.pcvar] = nipals(contents, dataPca.desiredComponents, dataPca.maxIterations);
      %         [dataPca.COEFF, dataPca.SCORE, dataPca.LATENT, dataPca.TSQUARED] = pca(contents, 'Algorithm', 'als');
      elapsed = toc;
      save(path, 'dataPca', 'columnLabels', 'elapsed');
      
    end
  end
  close all;
  fclose all;
  
end

% checking granger causality between individual datapoints is somewhat
% interesting.  everything granger causes the average
% value, which shouldn't be too much of a surpise.

%maybe something more
% interesting will come out of positive feedback granger causing eeg
% signals.  Or additional signals that aren't captured in this average
% value, like



