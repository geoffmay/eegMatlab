removeSham = true;
doPower = false;
removeSingleQuestions = true;
bootstrapRepetitions = 1;
debugMode = false;
rankOrder = true;

freqLabels = {'delta', 'theta', 'alpha', 'beta', 'hibeta'};
antLocs = antChannelLocs;
antLocs = antLocs(1:32);
if(~exist('sigTable', 'var') || debugMode)
  tic;
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
  if(removeSham)
    neuropsychData([6, 7], :) = [];
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
  numberOfEegMeasurements = (size(change.coherence, 1) * (size(change.coherence, 2) - 1)) / 2 * size(change.coherence, 3);
  effectSlopes = zeros(numberOfEegMeasurements, size(neuropsychData,2));
  pointCounts = zeros(numberOfEegMeasurements, size(neuropsychData,2));
  bootstrappedNeuropsych = cell(0);
  bootstrappedIndices = [];
  testableNeuropsychCounter = 0;
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
        
        subjectCount = size(neuropsychData, 1);
        testableNeuropsychCounter = testableNeuropsychCounter + 1;
        bootstrappedNeuropsych{testableNeuropsychCounter} = varName;
        bootstrappedIndices(testableNeuropsychCounter) = neuropsychColumn;
        if(~debugMode)
          for bootstrapCounter = 1:bootstrapRepetitions
            upperTriangleCounter = 0;
            fprintf('.');
            if(bootstrapCounter > 1)
              %ind = ceil((rand(subjectCount, 1)) .* subjectCount);
              ind = randperm(subjectCount);
              while(all(ind == 1:subjectCount))
                ind = randperm(subjectCount);
              end
            else
              ind = 1:subjectCount;
            end
            allInd(bootstrapCounter,:) = ind;
            y = y2 - y1;
            
            
            eegMeasureCounter = 1;
            isUpperTriangle = logical(zeros(1, numel(change.coherence)));
            for i1 = 1:size(change.coherence, 1)
              for i2 = 1:size(change.coherence, 2)
                if(i2 > i1)
                  upperTriangleCounter = upperTriangleCounter + 1;
                  isUpperTriangle(eegMeasureCounter) = 1;
                end
                for i3 = 1:size(change.coherence, 3)
                  if(bootstrapCounter == 1)
                    if(i2>i1)
                      allLabels{upperTriangleCounter} = sprintf('%s-%s %s', antLocs{i1}, antLocs{i2}, freqLabels{i3});
                    end
                  end
                  clear x yPerm;
                  yPerm = y(ind);
                  
                  for k = 1:size(neuropsychData, 1)
                    if(doPower)
                      x(k,1) = changes(k).abs(i,j);
                    else
                      x(k,1) = changes(k).coherence(i1,i2,i3);
                    end
                  end
                  remove = isnan(x) | isnan(y);
                  x(remove) = [];
                  yPerm(remove) = [];
                  yMax = mean(yPerm) + 1.5*std(yPerm);
                  yMin = mean(yPerm) - 1.5*std(yPerm);
                  xMax = mean(x) + 1.5*std(x);
                  xMin = mean(x) - 1.5*std(x);
                  
                  remove = (yPerm < yMin) | (yPerm > yMax) |  (x < xMin) | (x > xMax);
                  if(any(remove))
                    dummy = 1;
                    rhoMatrix(eegMeasureCounter, neuropsychColumn) = 0;
                    rhos(hitCounter,1) = 0;
                    ps(hitCounter,1) = 1;                    
                    npVarNames{hitCounter,1} = varName;
                    frequencies{hitCounter,1} = freqLabels{i3};
                    electrode1{hitCounter,1} = prePower.stat.channelLabels{i1};
                    electrode2{hitCounter,1} = prePower.stat.channelLabels{i2};
                    xs{hitCounter,1} = x;
                    ys{hitCounter,1} = yPerm;

                    hitCounter = hitCounter + 1;
                  else
                    if(length(x) > 2)
                      allCounter = allCounter + 1;
                      if(rankOrder)
                        [rho, p] = corr(x,yPerm);
                      else
                        [rho, p] = corr(x,yPerm, 'type', 'spearman');
                      end
                      pThreshold = 0.05;
                      pThreshold = 1.1;
                      if(p < pThreshold)
                        rhoMatrix(eegMeasureCounter, neuropsychColumn) = rho;
                        if(bootstrapCounter == 1)
                          rhos(hitCounter,1) = rho;
                          ps(hitCounter,1) = p;
                          npVarNames{hitCounter,1} = varName;
                          frequencies{hitCounter,1} = freqLabels{i3};
                          electrode1{hitCounter,1} = prePower.stat.channelLabels{i1};
                          electrode2{hitCounter,1} = prePower.stat.channelLabels{i2};
                          xs{hitCounter,1} = x;
                          ys{hitCounter,1} = yPerm;
                        end
                        hitCounter = hitCounter + 1;
                      end
                    end
                  end
                  eegMeasureCounter = eegMeasureCounter + 1;
                end
              end
            end
            %distribution testing
            stroopR = rhoMatrix(:, neuropsychColumn);
            flatStroop = stroopR(isUpperTriangle);
            fisherFlatStroop = FisherTransform(flatStroop);
            [thisT.H,thisT.P,thisT.CI,thisT.STATS] = ttest(fisherFlatStroop);
            ttestsNull(bootstrapCounter, testableNeuropsychCounter) = thisT;
            resampledMeanRs(bootstrapCounter, testableNeuropsychCounter) = mean(rhos);
            resampledStdRs(bootstrapCounter, testableNeuropsychCounter) = std(rhos);
            if(bootstrapCounter > 1)
              [thisT.H,thisT.P,thisT.CI,thisT.STATS] = ttest(fisherFlatStroop, originalFisherFlatStroop);
              ttestsResample(bootstrapCounter-1,testableNeuropsychCounter) = thisT;
            else
              originalFisherFlatStroop = fisherFlatStroop;
            end
          end
        end
      end
    end
  end
  
  
  
  sigTable = table(npVarNames, frequencies, electrode1, electrode2, rhos, ps, xs, ys);
  sigTable = sortrows(sigTable, 5);
  toc;
end

if(false)
  keep = logical(zeros(1, max(bootstrappedIndices)));
  keep(bootstrappedIndices) = 1;
  remove = ~keep;
  ttestsNull1 = ttestsNull;
  ttestsResample1 = ttestsResample;
  ttestsNull1(remove) = [];
  ttestsResample1(remove) = [];
end

if(bootstrapRepetitions > 1)
bootstrapSummary.ttestsNull = ttestsNull;
bootstrapSummary.ttestsResample = ttestsResample;
bootstrapSummary.resampledMeanRs = resampledMeanRs;
bootstrapSummary.resampledStdRs = resampledStdRs;
bootstrapSummary.effectSlopes = effectSlopes;
bootstrapSummary.measureLabels = allLabels;
bootstrapSummary.neuropsychLabels = bootstrappedNeuropsych;
bootstrapSummary.resampleHistory = allInd;
bootstrapSummary.note = 'the first iteration of each t-test is the original (non-permuted) version.  this is the first version that is permuted (resampled without replacement).';

save('/home/data/EEG/processed/bootstrapSummary2.mat', 'bootstrapSummary');

%interpretation of bootstrap results
resampleP = reshape([bootstrapSummary.ttestsResample.P], size(bootstrapSummary.ttestsResample));
resampleT = NaN(size(bootstrapSummary.ttestsResample));
for i = 1:size(bootstrapSummary.ttestsResample, 1)
  for j = 1:size(bootstrapSummary.ttestsResample, 2)
    resampleT(i,j) = bootstrapSummary.ttestsResample(i,j).STATS.tstat;
  end
end
greater = resampleT > 0;
lesser = resampleT < 0;
greaterSum = sum(greater, 1);
lesserSum = sum(lesser, 1);
imagesc(greater);
colorbar;
bootstrapSummary.neuropsychLabels(greaterSum > 90)

for i = 1:size(bootstrapSummary.ttestsNull, 2)
  for j = 2:size(bootstrapSummary.ttestsNull, 1)
    greaterCount(j-1, i) = bootstrapSummary.resampledMeanRs(1, i) > bootstrapSummary.resampledMeanRs(j, i);
  end
end
end

plotFilters = {'Multilingual_Aphasia_Adjusted', ...
  'STROOP_Int_Tscore', 'WAIS_Coding_Scaled', ...
  'WAIS_Digit_Scaled', 'WAIS_Symbol_Scaled',..., ...
  'Trail_Time_A',...,...
  'Trail_Errors_A',...,...
  'Trail_Time_B',...,...
  'Trail_Errors_B',...,...
  
  %neuropsych info
  %     'BSI_SOM_T_Score',...
  %     'BSI_OC_T_Score',...
  %     'BSI_IS_T_Score',...
  %     'BSI_DEP_T_Score',...
  %     'BSI_ANX_T_Score',...
  %     'BSI_HOS_T_Score',...
  %     'BSI_PHOB_T_Score',...
  %     'BSI_PAR_T_Score',...
  %     'BSI_PSY_T_Score',...
  %     'BSI_GSI_T_Score',...
  %     'BSI_PSDI_T_Score',...
  %     'BSI_PST_T_Score',...
  
  %%these are invididual items and need to be scored.
  %     'PCL5_a',...
  %     'PCL5_b',...
  %     'PCL5_c',...
  %     'PCL5_d',...
  %     'PCL5_d_additional1',...
  %     'PCL5_d_additional2',...
  %     'PCL5_d_additional3',...
  %     'PCL5_d_other',...
  %     'PCL5_e',...
  %     'PCL5_1',...
  %     'PCL5_2',...
  %     'PCL5_3',...
  %     'PCL5_4',...
  %     'PCL5_5',...
  %     'PCL5_6',...
  %     'PCL5_7',...
  %     'PCL5_8',...
  %     'PCL5_9',...
  %     'PCL5_10',...
  %     'PCL5_11',...
  %     'PCL5_12',...
  %     'PCL5_13',...
  %     'PCL5_14',...
  %     'PCL5_15',...
  %     'PCL5_16',...
  %     'PCL5_17',...
  %     'PCL5_18',...
  %     'PCL5_19',...
  %     'PCL5_20',...
  %     'PSQI_1',...
  %     'PSQI_2',...
  %     'PSQI_3',...
  %     'PSQI_4',...
  %     'PSQI_5_a',...
  %     'PSQI_5_b',...
  %     'PSQI_5_c',...
  %     'PSQI_5_d',...
  %     'PSQI_5_e',...
  %     'PSQI_5_f',...
  %     'PSQI_5_g',...
  %     'PSQI_5_h',...
  %     'PSQI_5_i',...
  %     'PSQI_5_j',...
  %     'PSQI_5_j_1',...
  %     'PSQI_6',...
  %     'PSQI_7',...
  %     'PSQI_8',...
  %     'PSQI_9',...
  %     'PSQI_10',...
  %     'PSQI_10_a',...
  %     'PSQI_10_b',...
  %     'PSQI_10_c',...
  %     'PSQI_10_d',...
  %     'PSQI_10_e',...
  %     'PSQI_10_e_1',...
  %     'Rivermead_a',...
  %     'Rivermead_b',...
  %     'Rivermead_c',...
  %     'Rivermead_d',...
  %     'Rivermead_e',...
  %     'Rivermead_f',...
  %     'Rivermead_g',...
  %     'Rivermead_h',...
  %     'Rivermead_i',...
  %     'Rivermead_j',...
  %     'Rivermead_k',...
  %     'Rivermead_l',...
  %     'Rivermead_m',...
  %     'Rivermead_n',...
  %     'Rivermead_o',...
  %     'Rivermead_p',...
  %     'Rivermead_q1',...
  %     'Rivermead_q1_answer',...
  %     'Rivermead_q2',...
  %     'Rivermead_q2_answer',...
  };





for plotFilterNumber = 1:length(plotFilters)
  try
    plotFilter = plotFilters{plotFilterNumber};
    if(length(plotFilter) > 0)
      stringInput = sigTable{:,1};
      match1 = strfind(sigTable{:,1},plotFilter);
      matches = find(cellfun(@length, match1));
    else
      matches = 1:size(sigTable,1);
    end
    
    doTopoplot = true;
    if(doTopoplot)
      newPs = sigTable{:,'ps'};
      significant = newPs(matches) < 0.05;
      plotThis = matches(significant);
      a = sigTable(plotThis,:);
      labels = cell(0);
      for i = 1:size(a,1)
        e1 = a{i,'electrode1'};
        e1 = e1{1};
        e2 = a{i,'electrode2'};
        e2 = e2{1};
        fr = a{i,'frequencies'};
        fr = fr{1};
        label = sprintf('%s-%s %s', e1, e2, fr);
        labels{i} = sprintf('%s-%s %s', e1, e2, fr);
      end
      labels = labels';
      plotVal = a{:,'rhos'};
      plotCoherencePca(plotVal, labels, 1, plotFilter);
      fprintf('\n%s (%d)', plotFilter, plotFilterNumber);
      %  export_fig(gcf,sprintf('%s.png',plotFilter),'-nocrop');
      
    end
  catch err
    %do nothing
    %rethrow(err);
  end
  close all;
end


if(false)
  
  bError = cellfun(@length, strfind(sigTable{:,1}, 'Trail_Errors_B')) > 0;
  theta = strcmp(sigTable{:,2}, 'theta');
  Fz = strcmp(sigTable{:,3}, 'Fz') | strcmp(sigTable{:,4}, 'Fz');
  Cz = strcmp(sigTable{:,3}, 'Cz') | strcmp(sigTable{:,4}, 'Cz');
  
  interesting = bError & theta & Fz & Cz;
  x = sigTable{interesting, 'xs'};
  y = sigTable{interesting, 'ys'};
  x = x{1};
  y = y{1};
  
  
end

if(false)
  
  
counter = 1;
clear noOutliers zMaxes;
for plotIndex = 1:length(matches)
  i = matches(plotIndex);
  close all;
  x = sigTable{i,7};
  y = sigTable{i,8};
  x = x{1};
  y = y{1};
  sty = std(y);
  mny = mean(y);
  clear z;
  for j = 1:length(x)
    z(j) = (y(j) - mny) / sty;
  end
  plotThis = true;
  if(max(abs(z)) < 1.5)
    noOutliers(counter) = i;
    zMaxes(counter) = max(abs(z));
    counter = counter + 1;
    plotThis = true;
    
  end
  if(plotThis)
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
end
end
%
% for counter = 1:length(noOutliers)
%   i = noOutliers(counter);
%   x = sigTable{i,7};
%   y = sigTable{i,8};
%   x = x{1};
%   y = y{1};
%
% end

