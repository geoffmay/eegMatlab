artifactOptions = {'Present', 'Removed'};

totalFileCounter = 1;

for artifactCounter = 1:length(artifactOptions)
  
  presentVsAbsent = artifactOptions{artifactCounter};
  
  %folder = '/Users/Geoff/Documents/reliability5/surface digest 2';
  %folder = '/Users/Geoff/Documents/reliability5/';
  %folder = '/Users/Geoff/Documents/reliability5/digest 3';
  %folder = '/Users/Geoff/Documents/MATLAB/processed/Oregon/reliability5';
  folder = ['/home/data/EEG/processed/Oregon/artifact', presentVsAbsent, 'Reliability'];
  %intermediateFolder = fullfile(folder, 'surface digest 2');
  intermediateFolder = ['/home/data/EEG/processed/Oregon/digest/artifact', presentVsAbsent, ''];
  %timeEstimateFolder = fullfile(folder, 'timeEstimates');
  timeEstimateFolder = ['/home/data/EEG/processed/Oregon/timeEstimates/artifact', presentVsAbsent, ''];
  
  mkdir(intermediateFolder);
  mkdir(timeEstimateFolder);
  
  files = dir(folder);
  files([files.isdir]) = [];
  files(cellfun(@length, strfind({files.name}, '.mat')) == 0) = [];
  
  doPlot = true;
  mat = NaN(70, length(files));
  fontSize= 20;
  
  %raw = load('~/Documents/reliability5/PM101.mat')
  %isCoh = cellfun(@length, strfind(raw.summary.surfLabels, 'coh')) > 0;
  
  %set fit options
  cohFitOptions = fitoptions('exp2');
  cohFitOptions.StartPoint = [0.0931   -0.0266    0.0948   -0.0009];
  cohFitOptions.Lower = [0   -1    0   -1];
  cohFitOptions.MaxIter = 2000;
  cohFitOptions.Upper = [3   0    3   0];
  powFitOptions = cohFitOptions;
  powFitOptions.Upper = [10   0    10   0];
  
  %create 'surface digest' by fitting curves
  stdCount = norminv(.95);
  for i = 1:length(files)
    fprintf('\n%d of %d', i, length(files));
    intermediateFile = fullfile(intermediateFolder, files(i).name);
    output.filenames(totalFileCounter) = intermediateFile;
    if(~exist(intermediateFile, 'file'))
      a = load(fullfile(folder, files(i).name));
      avgBlock = NaN(length(a.summary.surfResample), length(a.summary.surfLabels));
      stdBlock = NaN(length(a.summary.surfResample), length(a.summary.surfLabels));
      frameCounts = NaN(length(a.summary.surfResample),1);
      for windowCounter = 1:length(a.summary.surfResample)
        avgBlock(windowCounter,:) = a.summary.surfResample(windowCounter).averageDifference;
        stdBlock(windowCounter,:) = a.summary.surfResample(windowCounter).stddevDifference;
        frameCounts(windowCounter) = a.summary.surfResample(windowCounter).frameCount;
      end
      for measureCounter = 1:size(avgBlock,2)
        fprintf('\nmeasureCounter: %d', measureCounter);
        series95 = avgBlock(:, measureCounter) + (stdCount .* stdBlock(:, measureCounter));
        try
          [curveFit.FO, curveFit.G, curveFit.O] = fit(frameCounts, series95, 'exp2', cohFitOptions);
        catch ex
          curveFit.FO = ex;
          curveFit.G = ex;
          curveFit.O = ex;
        end
        fits(measureCounter) = curveFit;
      end
      save(intermediateFile, 'fits', '-v7.3');
    else
      load(intermediateFile);
    end
    
    
    
    if(false)
      %count max count for each coefficient
      clear coeffs;
      for i = 1:length(fits)
        try
          [coeffs(i, :)] = coeffvalues(fits(i).FO);
        catch ex
          [coeffs(i, :)] = NaN(1,4);
        end
      end
      for coeffIndex = 1:4
        tab = tabulate(coeffs(:, coeffIndex));
        [maxCount, index] = max(tab(:,2));
        commonCoefficientsCounts(totalFileCounter, coeffIndex) = maxCount;
        commonCoefficients(totalFileCounter, coeffIndex) =  tab(index, 1);
      end
      totalFileCounter = totalFileCounter + 1;
    end
  end
end


output.commonCoefficientsCounts = commonCoefficientsCounts;
output.commonCoefficients = commonCoefficients; 
output.notes = 'columns are a,b,c,d in the equation y = a*x^b + c*x^d';

save('/home/data/EEG/processed/Oregon/phaseFitCounts', 'output', '-v7.3');

channel = 0;
x = (((1:size(topography.estimatedTimeLag, 1)) - 1024 )./ 1024)';

channel = channel + 1;
close all;
for i = 1:size(topography.frequencyLimits, 1);
  figure;
  plot(x, topography.estimatedTimeLag(:, channel, i));
  title(sprintf('%s %d-%d', topography.chanlocs(channel).labels, topography.frequencyLimits(i, 1), topography.frequencyLimits(i, 2)));
end
tilefigs

close all;
x = (((1:size(topography.estimatedTimeLag, 1)) - 1024 )./ 1024)';
for i = 1:size(topography.estimatedTimeLag, 2);
  figure;
  plot(x, topography.estimatedTimeLag(:, i, 1));
  title(sprintf('%s %d-%d', topography.chanlocs(i).labels, topography.frequencyLimits(1, 1), topography.frequencyLimits(1, 2)));
  ylim([-0.005, 0.0055]);
end
tilefigs

%   norms = load('/Users/Geoff/Documents/MATLAB/processed/Oregon/avgValuesProcessed/subjectNorms.mat');
%   syms x;
%
%
%   %create time estimates
%   stdFractions = [.001 .002 .005 .01 .02 .05 .1 .2 .5 1];
%   data = load(fullfile(intermediateFolder, files(1).name));
%
%   for i = 1:length(files)
%     timeEstimateFile = fullfile(timeEstimateFolder, files(i).name);
%     if(~exist(timeEstimateFile))
%       data = load(fullfile(intermediateFolder, files(i).name));
%       for j = 1:length(data.fits)
%         fprintf('\nfile %d of %d, measure %d of %d', i, length(files), j, length(data.fits));
%         coeffs = coeffvalues(data.fits(j).fun);
%         fasterFirst = abs(coeffs(2)) > abs(coeffs(4));
%         if(~fasterFirst)
%           temp = coeffs(1:2);
%           coeffs(1:2) = coeffs(3:4);
%           coeffs(3:4) = temp;
%         end
%         coeffMat(i, j, :) = coeffs;
%         a = coeffs(1);
%         b = coeffs(2);
%         c = coeffs(3);
%         d = coeffs(4);
%         for timeEstimateCounter = 1:length(stdFractions)
%           clear timeEstimate;
%           stdFraction = stdFractions(timeEstimateCounter);
%           %timeEstimate.stdFraction = stdFractions(timeEstimateCounter);
%           timeEstimate.criterion = sprintf('time to reach %f * std dev(population mean)', stdFractions(timeEstimateCounter));
%           %                 timeEstimate.durations = NaN(length(files), length(data.fits));
%           y = norms.summary.intersubjectStd(j) * stdFraction;
%           eqn = a*exp(x*b) + c*exp(x*d) == y;
%           sol = solve(eqn);
%           %                 timeEstimate.durations(i,j) = double(sol);
%           timeEstimate.durations = double(sol);
%           timeEstimate.criterion = sprintf('time to reach difference < %f * std dev(population mean)', stdFractions(timeEstimateCounter));
%           timeEstimates(j, timeEstimateCounter, 1) = timeEstimate;
%
%           y = stdFraction;
%           eqn = a*exp(x*b) + c*exp(x*d) == y;
%           sol = solve(eqn);
%           %                 timeEstimate.durations(i,j) = double(sol);
%           timeEstimate.durations = double(sol);
%           timeEstimate.criterion = sprintf('time to reach difference < %f)', stdFractions(timeEstimateCounter));
%           timeEstimates(j, timeEstimateCounter, 2) = timeEstimate;
%         end
%       end
%       measureLabels = norms.summary.labels;
%       save(timeEstimateFile, 'timeEstimates', 'measureLabels');
%     end
%   end
%
%   %compute mean estimate
%   highPrecision = false;
%   timeEstimateFiles = dir(timeEstimateFolder);
%   timeEstimateFiles([timeEstimateFiles.isdir]) = [];
%   timeEstimateCount = 0;
%   figFolder = '/Users/Geoff/Documents/figures/reliability';
%   for i = 1:length(timeEstimateFiles)
%     fprintf('\nfile %d of %d', i, length(timeEstimateFiles));
%     data = load(fullfile(timeEstimateFolder, timeEstimateFiles(i).name));
%     if(isfield(data, 'timeEstimates'))
%       timeEstimates = data.timeEstimates;
%       if(~highPrecision)
%         tes = [timeEstimates(:,9,1).durations]; %9th column = 0.5 std devs
%       else
%         tes = [timeEstimates(:,6,1).durations]; %9th column = 0.05 std devs
%       end
%       [tes1.data, tes1.labels] = pruneCoherenceMeasures(tes, measureLabels);
%
%       timeEstimates = log([tes1.data]);
%       timeEstimateMatrix(i, :) = [tes1.data];
%       %timeEstimateMatrix(i, :) = [tes1.data];
%       if(timeEstimateCount == 0)
%         timeEstimateSum = timeEstimates;
%       else
%         timeEstimateSum = timeEstimateSum + timeEstimates;
%       end
%       timeEstimateCount = timeEstimateCount + 1;
%     end
%   end
%   timeEstimateMeanLog = timeEstimateSum ./ timeEstimateCount;
%   timeEstimateMedian = median(timeEstimateMatrix);
%   timeEstimateMean = exp(timeEstimateMeanLog);
%   %trimmed = timeEstimateMean(timeEstimateMean < 3000);  hist(trimmed);
%   figure;
%
%   if(highPrecision)
%     figTitle = 'time to reach precision of 0.05 population standard deviations';
%     trimmed = timeEstimateMedian(timeEstimateMedian < (150*60));
%     hist(trimmed./60);
%     title(figTitle)
%     xlabel('median time to reach convergence(minutes)');
%     ylabel('power and coherence measure count');
%     set(gcf,'color','white');
%     export_fig(fullfile(figFolder, 'hist.png'));
%     plotCoherencePca(timeEstimateMedian ./ 60, tes1.labels, [40 120], figFolder);
%   else
%     figTitle = 'time to reach precision of 0.5 population standard deviations';
%     trimmed = timeEstimateMedian(timeEstimateMedian < (60*60));
%     hist(trimmed./60);
%     title(figTitle)
%     xlabel('median time to reach convergence(minutes)');
%     ylabel('power and coherence measure count');
%     set(gcf,'color','white');
%     export_fig(fullfile(figFolder, 'hist.png'));
%     plotCoherencePca(timeEstimateMedian ./ 60, tes1.labels, [0 45], figFolder);
%   end
%
%
%   if(false)
%     %demonstrate
%     label = 'coh Fp1-Fp2 5Hz-8Hz';
%     ind = find(strcmp(measureLabels, label));
%     te = squeeze(timeEstimates(ind, :, :));
%     tes = [timeEstimates(:,7,1).durations];
%     [tes1.data, tes1.labels] = pruneCoherenceMeasures(tes, measureLabels);
%     %plotCoherencePca(log(tes1.data), tes1.labels);
%     dat = tes1.data ./ 60;
%     close all;
%     plotCoherencePca(dat, tes1.labels, [5 120]);
%
%     isPower = cellfun(@length, strfind(tes1.labels, 'abs')) > 0;
%     isTheta = cellfun(@length, strfind(tes1.labels, '5Hz-8Hz')) > 0;
%     ind = find(isPower & isTheta);
%     thetaPower = tes1.data(ind)
%     for i = 1:length(thetaPower)
%       fprintf('\n%s: %0.0f minutes', tes1.labels{ind(i)}, thetaPower(i) / 60);
%     end
%   end
%
%   if(false)
%     %save('exponentialFit.mat', 'coeffMat');
%     % equation is a*e^(b*x) + c*e^(d*x)
%     %   'coh icaComp1-icaComp19 9Hz-12Hz' is frickin enormous for a and c
%     %   (9e12 and -9e12)
%     maxVal = max(max(max(coeffMat)));
%     minVal = min(min(min(coeffMat)));
%     maxInd = find(coeffMat == maxVal);
%     minInd = find(coeffMat == minVal);
%     minVal = min(min(min(coeffMat)));
%
%     extreme = abs(coeffMat) > .5;
%     coeffMat(extreme)  = 0;
%     cohMatA = squeeze(coeffMat(:,isCoh,1));
%     powMatA = squeeze(coeffMat(:,~isCoh,1));
%     cohMatB = squeeze(coeffMat(:,isCoh,2));
%     powMatB = squeeze(coeffMat(:,~isCoh,2));
%     cohMatC = squeeze(coeffMat(:,isCoh,3));
%     powMatC = squeeze(coeffMat(:,~isCoh,3));
%     cohMatD = squeeze(coeffMat(:,isCoh,4));
%     powMatD = squeeze(coeffMat(:,~isCoh,4));
%
%     imagesc(cohMatB);
%     colorbar;
%
%     meanD = mean(cohMatD, 2);
%     [~, indD] = sort(meanD);
%     sortCohMatD = cohMatD(indD,:);
%     imagesc(sortCohMatD);
%
%     if(false)
%       %find extremes
%       %here is a very low d (low non-stationarity?)
%       minD = cohMatD == min(min(cohMatD));
%       minD1 = find(any(minD,2));
%       minD2 = find(any(minD,1));
%       i = minD1;
%       ind = minD2;
%       a = load(fullfile(folder, files(i).name));
%       b = a.plot95Raw{ind};
%
%       %here is a very low b (low measurement noise?)
%       minB = cohMatB == min(min(cohMatB));
%
%       minB1 = find(any(minB,2));
%       minB2 = find(any(minB,1));
%       i = minB1;
%       ind = minB2;
%       a = load(fullfile(folder, files(i).name));
%       b = a.plot95Raw{ind};
%
%
%     end
%
%     %check for neuropsych correlations
%     neuropsychDataFilename = '~/Documents/MATLAB/EEG/PTSD MIND for Geoff.mat';
%     neuropsychData = load(neuropsychDataFilename);
%     neuropsychData = neuropsychData.neuropsychData;
%
%     for i = 1:length(files)
%       [folder, file, ext] = fileparts(files(i).name);
%       file = str2num(file(4:5));
%       rows(i) = find(neuropsychData{:,1} == file);
%     end
%     for i = 1:size(neuropsychData, 2)
%       x = meanD;
%       y = neuropsychData{rows, i};
%       [rhos(i), ps(i)] = corr(x,y);
%     end
%     pInd = ps < 0.05
%
%
%     hist(cohMat);
%     legend({'a', 'b', 'c', 'd'});
%   end
%
% end
%

