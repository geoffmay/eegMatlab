doPower = false;
removeSingleQuestions = true;

if(false)
%load neuropsychdata
neuropsychFilename = '/home/data/EEG/processed/Robi/robiNeuropsych.mat';
if(~exist(neuropsychFilename, 'file'))
  neuropsychData = readcsv('/home/data/EEG/processed/Robi/ROBI Database_Pre_Post.csv');
%   pre = readcsv('/home/data/EEG/processed/Robi/ROBI Neuropsych_Pre.csv');
%   post = readcsv('/home/data/EEG/processed/Robi/ROBI Neuropsych_Post.csv');
  save(fullfile(folder,'robiNeuropsych.mat'),'neuropsychData');
end
load(neuropsychFilename);
if(removeSingleQuestions)
  varNames = neuropsychData.Properties.VariableNames;
  remove = zeros(1,length(varNames));
  filter = [{'Rivermead'},{'PSQI'},{'PCL5'},{'BSI_1'},{'BSI_2'},{'BSI_3'},...
    {'BSI_4'},{'BSI_5'},{'BSI_6'},{'BSI_7'},{'BSI_8'},{'BSI_9'}];
  for i = 1:length(filter)
    hits = cellfun(@length, strfind(varNames, filter{i}));
    remove = remove | hits;
  end
  neuropsychData(:, remove) = [];
end

freqLabels = [{'delta'},{'theta'},{'alpha'},{'beta'},{'hibeta'}];

%load EEG power into memory.
if(doPower)
  powerFolder = '/home/data/EEG/processed/Robi/power';
  powerFiles = dir(powerFolder);
  powerFiles = {powerFiles.name};
else
  powerFolder = '/home/data/EEG/processed/Robi/coherenceReref';
  powerFiles = dir(powerFolder);
  powerFiles = {powerFiles.name};
end

for subjectNumber = 1:size(neuropsychData,1)
  subjectId = neuropsychData{subjectNumber, 'Subj_ID'};
  subjectId = subjectId{1};
  idNumber = subjectId(5:7);
  searchTarget1 = sprintf('ROBI_%s_baseline eyes open', idNumber);
  searchTarget2 = sprintf('ROBI_%s_outcome eyes open', idNumber);
  preFilenumber = find(cellfun(@length, strfind(powerFiles, searchTarget1)));
  postFilenumber = find(cellfun(@length, strfind(powerFiles, searchTarget2)));
  prePower = load(fullfile(powerFolder, powerFiles{preFilenumber}));
  postPower = load(fullfile(powerFolder, powerFiles{postFilenumber}));
  if(doPower)
    change.abs = -prePower.stat.meanAbsPowers + postPower.stat.meanAbsPowers;
    change.rel = -prePower.stat.meanRelPowers + postPower.stat.meanRelPowers;
    change.asym = -prePower.stat.meanAsymmetries + postPower.stat.meanAsymmetries;
    change.ratios = -prePower.stat.meanPowerRatios + postPower.stat.meanPowerRatios;
  else
    change.coherence =postPower.stat.meanCoherences - prePower.stat.meanCoherences;
    change.coherenceStd =postPower.stat.stdDevCoherences - prePower.stat.stdDevCoherences;
  end
  changes(subjectNumber) = change;
end



hitCounter = 1;
allCounter = 0;
for neuropsychColumn = 1:size(neuropsychData,2)
  y1 = neuropsychData{:,neuropsychColumn};
  if(isnumeric(y1))
    varName = neuropsychData.Properties.VariableNames{neuropsychColumn};
    searchName = sprintf('%s_Post', varName);
    postIndex = find(strcmp(neuropsychData.Properties.VariableNames, searchName));
    if(length(postIndex) > 0)
      fprintf('\n%s',varName);
      y2 = neuropsychData{:,postIndex(1)};
      if(~isnumeric(y2))
        for(rewrite = 1:length(y2))
          y2N = str2num(y2{rewrite});
          if(length(y2N) == 0)
            y2N = NaN;
          end
          y2Num(rewrite,1) = y2N;
        end
        y2 = y2Num;
      end
      for i1 = 1:size(change.coherence, 1)
        for i2 = 1:size(change.coherence, 2)
          for i3 = 1:size(change.coherence, 3)
            clear x y;
            y = y2 - y1;
            for k = 1:size(neuropsychData, 1)
              if(doPower)
                x(k,1) = changes(k).abs(i,j);
              else
                x(k,1) = changes(k).coherence(i1,i2,i3);
              end
            end
            remove = isnan(x) | isnan(y);
            x(remove) = [];
            y(remove) = [];
            yMax = mean(y) + 3*std(y);
            yMin = mean(y) - 3*std(y);
            xMax = mean(x) + 3*std(x);
            xMin = mean(x) - 3*std(x);
            
            remove = (y < yMin) | (y > yMax) |  (x < xMin) | (x > xMax);
            if(any(remove))
              dummy = 1;
            end
            
            
            if(length(x) > 2)
              allCounter = allCounter + 1;
              [rho, p] = corr(x,y);
              if(p < 0.05)
                %               close all;
                %               scatter(x,y);
                rhos(hitCounter,1) = rho;
                ps(hitCounter,1) = p;
                varNames{hitCounter,1} = varName;
                frequencies{hitCounter,1} = freqLabels{i3};
                electrode1{hitCounter,1} = prePower.stat.channelLabels{i1};
                electrode2{hitCounter,1} = prePower.stat.channelLabels{i2};
                xs{hitCounter,1} = x;
                ys{hitCounter,1} = y;
                hitCounter = hitCounter + 1;
              end
            end
          end
        end
      end
    end
  end
end

sigTable = table(varNames, frequencies, electrode1, electrode2, rhos, ps, xs, ys);
sigTable = sortrows(sigTable, 5);
end


plotFilter = 'WAIS';
if(length(plotFilter) > 0)
  match1 = strfind(sigTable{:,1},plotFilter);
  matches = find(cellfun(@length, match1));
else
  matches = 1:size(sigTable,1);
end

for plotIndex = 1:length(matches)
  i = matches(plotIndex);
  close all;
  x = sigTable{i,7};
  y = sigTable{i,8};
  x = x{1};
  y = y{1}; 
  figure;
  scatter(x, y);
  neuropsychMeasure = sigTable{i,1};
  freq = sigTable{i,2};
  electrd1 = sigTable{i,3};
  electrd2 = sigTable{i,4};
  neuropsychMeasure = neuropsychMeasure{1};
  neuropsychMeasure = strrep(neuropsychMeasure, '_', ' ');
  freq = freq{1};
  electrd1 = electrd1{1};
  electrd2 = electrd2{1};
  pValue = sigTable{i,6};
  rValue = sigTable{i,5};
  title(sprintf('R=%d p=%d',rValue,pValue));
  xlabel(sprintf('%s %s-%s coherence change', freq, electrd1, electrd2));
  ylabel(sprintf('change in %s', neuropsychMeasure));
end

