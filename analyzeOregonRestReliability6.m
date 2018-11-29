%folder = '/Users/Geoff/Documents/reliability5/surface digest 2';
solveForY = false;

removeEyeblink = false;
printIndex = 1;
targetStd = [0.05];
figCounter = 101;    inputFolders = {'/home/data/EEG/processed/Oregon/artifactPresentReliability', '/home/data/EEG/processed/Oregon/artifactRemovedReliability'};
    inputFileCounter = 1;
    outputFolder = '/home/data/EEG/processed/Oregon/reliabilityFits';
    close all;
    if(~exist(workingFile, 'file'))
      for i = 1:length(inputFolders)
        files = dir(inputFolders{i});
        inputFiles(inputFileCounter:(inputFileCounter + length(files) - 1)) = files;
        correspondingInputFolders(inputFileCounter:(inputFileCounter + length(files) - 1)) = inputFolders(i);
        inputFileCounter = inputFileCounter + length(files);
      end
      remove = false(1, length(inputFiles));
      % inputFiles([inputFiles.isdir]) = [];
      % inputFiles(cellfun(@length, strfind({inputFiles.name}, '.mat')) == 0) = [];
      remove([inputFiles.isdir]) = 1;
      remove(cellfun(@length, strfind({inputFiles.name}, '.mat')) == 0) = 1;
      inputFiles(remove) = [];
      correspondingInputFolders(remove) = [];


idByGroup = [1, 2;2, 2;3, 1;4, 3;5, 3;6, 3;7, 3;8, 3;9, 2;10, 3;11, 3;12, 3;15, 2;16, 3;17, 1;18, 1;19, 1;20, 1;21, 1;22, 2;23, 2;24, 1;25, 1;26, 1;27, 1;28, 1;29, 1;30, 1;31, 1;32, 3;33, 3;34, 3;35, 3;36, 3;37, 1;38, 2;39, 2;41, 3;42, 2;43, 2;44, 2;45, 2;46, 2;47, 2;48, 2];


for targetStd = [0.05, 0.1, 0.5]
  for removeEyeblink = [true, false]
    
    imageType = '-dpng';
    workingFile = '/home/data/EEG/processed/Oregon/projectedZScoresBulk.mat';
    inputFolders = {'/home/data/EEG/processed/Oregon/artifactPresentReliability', '/home/data/EEG/processed/Oregon/artifactRemovedReliability'};
    inputFileCounter = 1;
    outputFolder = '/home/data/EEG/processed/Oregon/reliabilityFits';
    close all;
    if(~exist(workingFile, 'file'))
      for i = 1:length(inputFolders)
        files = dir(inputFolders{i});
        inputFiles(inputFileCounter:(inputFileCounter + length(files) - 1)) = files;
        correspondingInputFolders(inputFileCounter:(inputFileCounter + length(files) - 1)) = inputFolders(i);
        inputFileCounter = inputFileCounter + length(files);
      end
      remove = false(1, length(inputFiles));
      % inputFiles([inputFiles.isdir]) = [];
      % inputFiles(cellfun(@length, strfind({inputFiles.name}, '.mat')) == 0) = [];
      remove([inputFiles.isdir]) = 1;
      remove(cellfun(@length, strfind({inputFiles.name}, '.mat')) == 0) = 1;
      inputFiles(remove) = [];
      correspondingInputFolders(remove) = [];
      
      doPlot = true;
      mat = NaN(70, length(inputFiles));
      fontSize= 20;
      
      raw = load('/home/data/EEG/processed/Oregon/artifactPresentReliability/PM101Surface.mat');
      isCoh = cellfun(@length, strfind(raw.summary.surfLabels, 'coh')) > 0;
      isAsy = cellfun(@length, strfind(raw.summary.surfLabels, 'asy')) > 0;
      isAbs = cellfun(@length, strfind(raw.summary.surfLabels, 'abs')) > 0;
      
      
      %perform double-exponential fitting
      linearFitOptions = fitoptions('poly1');
      
%Example fit type      
%       ft = fittype( @(a, b, c, d, x, y) a*x.^2+b*x+c*exp(...
%         -(y-d).^2 ), 'independent', {'x', 'y'},...
%         'dependent', 'z' );

      
      inverseFitOptions = fitoptions('poly1');
      exp1FitOptions = fitoptions('exp1');
      
      cohFitOptions = fitoptions('exp2');
      cohFitOptions.StartPoint = [0.0931   -0.0266    0.0948   -0.0009];
      cohFitOptions.Lower = [0   -1    0   -1];
      cohFitOptions.MaxIter = 2000;
      cohFitOptions.Upper = [3   0    3   0];
      
      powFitOptions = fitoptions('exp2');
      powFitOptions.StartPoint = [0.0931   -0.0266    0.0948   -0.0009];
      powFitOptions.Lower = [0   -1    0   -1];
      powFitOptions.MaxIter = 2000;
      powFitOptions.Upper = [10   0    10   0];
      asyFitOptions = fitoptions('exp2');
      asyFitOptions.StartPoint = [0.0931   -0.0266    0.0948   -0.0009];
      asyFitOptions.Lower = [0   -1    0   -1];
      asyFitOptions.MaxIter = 2000;
      asyFitOptions.Upper = [1000   0    1000   0];
      stdCoeff = norminv(.95);
      for i = 1:length(inputFiles)
        [partFolder, partFile, partExt] = fileparts(correspondingInputFolders{i});
        outfile = fullfile(outputFolder, partFile, inputFiles(i).name);
        if(~exist(outfile, 'file'))
          placeholder = sprintf('started on %s', char(datetime));
          save(outfile, 'placeholder');
          a = load(fullfile(correspondingInputFolders{i}, inputFiles(i).name));
          searchString = 'coh Fp1-Fp2 5Hz-8Hz';
          plotInd = find(strcmp(a.summary.surfLabels, searchString));
          
          %    ind = find(strcmp({a.summaries.label}, 'coh Fp1-Fp2 5Hz-8Hz'));
          ind = find(strcmp(a.summary.surfLabels, 'coh Fp1-Fp2 5Hz-8Hz'));
          %    ind = 1;
          row = mat(:, i);
          plotInd = find(strcmp(a.summary.surfLabels, searchString));
          avgMatrix = NaN(length(a.summary.surfResample), length(a.summary.surfResample(1).averageDifference));
          conf95Matrix = NaN(size(avgMatrix));
          conf95MatrixAlt = NaN(size(avgMatrix));
          for j = 1:length(a.summary.surfResample)
            thisAvg = a.summary.surfResample(j).averageDifference;
            thisStd = a.summary.surfResample(j).stddevDifference;
            avgMatrix(j, :) = thisAvg;
            conf95Matrix(j, :) = thisAvg + (stdCoeff .* a.summary.surfResample(j).stddevDifference);
            ind95 = find(a.summary.surfResample(j).percentiles.percentileKeys == 95);
            conf95MatrixAlt(j, :) = a.summary.surfResample(j).percentiles.percentileValues(ind95,:);
          end
          
          b = conf95MatrixAlt(:, ind)';
          x = (1:length(b)) .* 15;
          %     row(1:length(b)) = b;
          %     mat(:,i) = row;
          
          clear fits
          for ind = 1:size(conf95MatrixAlt, 2)
            fprintf('\nfile %d of %d; measure %d of %d', i, length(inputFiles), ind, size(conf95MatrixAlt, 2));
            b = conf95MatrixAlt(:, ind)';
            x = (1:length(b)) .* 15;
            %         row(1:length(b)) = b;
            %         mat(:,i) = row;
            if(~any(isnan(b)))
              if(isCoh(ind))
                [f.fun, f.good] = fit(x', b', 'exp2', cohFitOptions);
              elseif(isAsy(ind))
                [f.fun, f.good] = fit(x', b', 'exp2', asyFitOptions);
              else
                [f.fun, f.good] = fit(x', b', 'exp2', powFitOptions);
              end
              if(ind == 1)for targetStd = [0.05, 0.1, 0.5]
                  for removeEyeblink = [true, false]
                    
                    fits = repmat(f, [1, size(conf95MatrixAlt, 2)]);
                    goods = NaN(1, size(conf95MatrixAlt, 2));
                  end
                  fits(ind) = f;
                  goods(ind) = f.good.adjrsquare;
                  if(doPlot)
                    if(ind == plotInd)
                      close all;
                      figure;
                      scatter(x, b);
                      hold on;
                      plot(f.fun);
                      label = a.summary.surfLabels{ind};
                      label = 'Difference between subsamples of theta power at Fp1';
                      title(label, 'fontsize', fontSize);
                      %             legend('95 %ile');
                      legend({'95 %ile', sprintf('fit (R^2 = %f)', f.good.adjrsquare)});
                      
                      ylabelCoherence = 'difference between sample means (% coherence)';
                      ylabelPower = 'difference between sample means (log(\muV^2))';
                      ylabel(ylabelPower, 'fontsize', fontSize);
                      xlabel('sample duration (seconds)', 'fontsize', fontSize);
                      % title('Difference between subsamples of theta coherence at Fp1-Fp2', 'fontsize', fontSize);
                      set(gca, 'fontsize', fontSize);
                      set(gcf, 'color', 'white');
                    end
                  end
                end
              end
              labels = a.summary.surfLabels;
              save(outfile, 'fits', 'goods', 'labels', '-v7.3');
            end
          end
          
          %example variation over time
          if(false)
            load('/home/data/EEG/processed/Oregon/artifactPresent/PM101.mat');
            a = surCoh.matrix(:, find(strcmp(surCoh.labels, 'abs Oz 9Hz-12Hz')));
            plot((1:length(a)) ./ 128, a);
            title('Oz alpha power over time');
            xlabel('sample duration (seconds)')
            ylabel('alpha power (log(\muV^2))');
          end
          
          %show R values
          rejectArtifact = false;
          if(rejectArtifact)
            keepFiles = strcmp(correspondingInputFolders, '/home/data/EEG/processed/Oregon/artifactRemovedReliability');
          else
            keepFiles = ~strcmp(correspondingInputFolders, '/home/data/EEG/processed/Oregon/artifactRemovedReliability');
          end
          filterLabels = {'EXG', 'Ana'}; %exclude these
          if(~exist('good', 'var'))
            load(outfile);
          end
          keepLabelInd = true(1, length(goods));
          a = load(fullfile(correspondingInputFolders{1}, inputFiles(1).name));
          for i = 1:length(filterLabels)
            hit = cellfun(@length, strfind(a.summary.surfLabels, filterLabels{i})) > 0;
            keepLabelInd(hit) = 0;
          end
          for i = 1:length(keepLabelInd)
            parts = strsplit(a.summary.surfLabels{i}, ' ');
            if(strcmp(parts{1}, 'asy'))
              locs = strsplit(parts{2}, '-');
              if(strcmp(locs{1}, locs{2}))
                keepLabelInd(i) = 0;
              end
            end
          end
          keepLabels = a.summary.surfLabels(keepLabelInd);
          clear adjRSquaredMatrix;
          rSquaredCounter = 1;
          inputFiles = inputFiles(keepFiles);
          correspondingInputFolders = correspondingInputFolders(keepFiles);
          for i = 1:length(inputFiles)
            fprintf('\nfile %d of %d', i, length(inputFiles));
            [partFolder, partFile, partExt] = fileparts(correspondingInputFolders{i});
            outfile = fullfile(outputFolder, partFile, inputFiles(i).name);
            if(exist(outfile, 'file'))
              fileContents = load(outfile);
              goodRow = fileContents.goods(keepLabelInd);
              keepInd = ~isnan(goodRow);
              keepLabels1 = keepLabels(keepInd);
              adjRSquaredMatrix(rSquaredCounter,:) = goodRow(keepInd);
              rSquaredCounter = rSquaredCounter + 1;
            end
            % figure;
            % hist(goods(keepLabel));
          end
          % hist(reshape(adjRSquaredMatrix, [1 numel(adjRSquaredMatrix)]), 100);
          % title('Goodness of fit for y = ax^b + cx^d');
          % xlabel('adjusted R squared');
          % ylabel('number of qEEG variables');
          
          
          valuesFilename = '/home/data/EEG/processed/Oregon/ptsdValues.mat';
          
          if(solveForY)
            projectedZScores = '/home/data/EEG/processed/Oregon/projectedZScores.mat';
          else
            projectedZScores = '/home/data/EEG/processed/Oregon/projectedZScoresBulk.mat';
          end
          
          
          
          if(~exist(valuesFilename, 'file'))
            clear norms;
            valuesFolders = {'/home/data/EEG/processed/Oregon/artifactPresentValues', '/home/data/EEG/processed/Oregon/artifactRemovedValues'};
            for folderCounter = 1:length(valuesFolders)
              valuesFiles = dir(valuesFolders{folderCounter});
              valuesFiles([valuesFiles.isdir]) = [];
              for i = 1:length(valuesFiles)
                fprintf('\n%s, %d of %d, %s', char(datetime), i, length(valuesFiles), valuesFolders{folderCounter});
                a = load(fullfile(valuesFolders{folderCounter}, valuesFiles(i).name));
                %       if(strcmp(valuesFolders{folderCounter}, '/home/data/EEG/processed/Oregon/artifactPresentValues'))
                %         norms.artifactPresent.filenames(i) = fullfile(valuesFolders{folderCounter}, valuesFiles(i).name);
                %       else
                %         norms.artifactRemoved.filenames(i) = fullfile(valuesFolders{folderCounter}, valuesFiles(i).name);
                for j = 1:length(a.summaries)
                  if(i == 1)
                    if(strcmp(valuesFolders{folderCounter}, '/home/data/EEG/processed/Oregon/artifactPresentValues'))
                      norms.artifactPresentLabels{j} = a.summaries(j).measureLabel;            %       end
                      
                    else
                      norms.artifactRemovedLabels{j} = a.summaries(j).measureLabel;
                    end
                  end
                  if(strcmp(valuesFolders{folderCounter}, '/home/data/EEG/processed/Oregon/artifactPresentValues'))
                    norms.artifactPresentMeans(i,j) = a.summaries(j).meanValue;
                    norms.artifactPresentStddevs(i,j) = a.summaries(j).stdValue;
                  else
                    norms.artifactRemovedMeans(i,j) = a.summaries(j).meanValue;
                    norms.artifactRemovedStddevs(i,j) = a.summaries(j).stdValue;
                  end
                end
              end
            end
            
            %restructure data
            temp = norms;
            clear norms;
            norms.rawData = temp;
            
            %track source filenames
            for folderCounter = 1:length(valuesFolders)
              valuesFiles = dir(valuesFolders{folderCounter});
              valuesFiles([valuesFiles.isdir]) = [];
              for i = 1:length(valuesFiles)
                fprintf('\n%s, %d of %d, %s', char(datetime), i, length(valuesFiles), valuesFolders{folderCounter});
                if(strcmp(valuesFolders{folderCounter}, '/home/data/EEG/processed/Oregon/artifactPresentValues'))
                  norms.rawData.artifactPresentFilenames{i} = fullfile(valuesFolders{folderCounter}, valuesFiles(i).name);
                else
                  norms.rawData.artifactRemovedFilenames{i} = fullfile(valuesFolders{folderCounter}, valuesFiles(i).name);
                end
              end
            end
            
            for i = 1:size(norms.artifactPresentMeans, 2)
              artifactPresent.means = mean(norms.rawData.artifactPresentMeans);
              artifactPresent.intersubjectStddev = std(norms.rawData.artifactPresentMeans);
              artifactPresent.intrasubjectStddev = mean(norms.rawData.artifactPresentStddevs);
              artifactPresent.measureLabels = norms.rawData.artifactPresentLabels;
              artifactPresent.means = mean(norms.rawData.artifactPresentMeans);
              
              artifactRemoved.means = mean(norms.rawData.artifactRemovedMeans);
              artifactRemoved.intersubjectStddev = std(norms.rawData.artifactRemovedMeans);
              artifactRemoved.intrasubjectStddev = mean(norms.rawData.artifactRemovedStddevs);
              artifactRemoved.measureLabels = norms.rawData.artifactRemovedLabels;
              artifactRemoved.means = mean(norms.rawData.artifactRemovedMeans);
              
              norms.artifactPresent = artifactPresent;
              norms.artifactRemoved = artifactRemoved;
            end
            save(valuesFilename, 'norms', '-v7.3');
          else
            norms = load(valuesFilename);
            norms = norms.norms;
          end
          
          %fit a number of Z scores.
          if(solveForY)
            %   zSamples = 0.01 : 0.01 : 2;
            zSamples = [0.5, 0.05];
            s.zSamples = zSamples;
            kCount = length(zSamples);
          else
            times = 15:15:(60 * 120); %assume 120 minutes is enough?
            s.times = times;
            kCount = length(times);
          end
          fitFolders = {'/home/data/EEG/processed/Oregon/reliabilityFits/artifactRemovedReliability', '/home/data/EEG/processed/Oregon/reliabilityFits/artifactPresentReliability'};
          for folderCounter = 1:length(fitFolders)
            fitFiles = dir(fitFolders{folderCounter});
            fitFiles([fitFiles.isdir]) = [];
            for i = 1:length(fitFiles)
              fit = load(fullfile(fitFolders{folderCounter}, fitFiles(i).name));
              if(strcmp(fitFolders{folderCounter}, '/home/data/EEG/processed/Oregon/reliabilityFits/artifactRemovedReliability'))
                thisNorm = norms.artifactRemoved;
                if(exist('projections', 'var'))
                  if(isfield(projections, 'artifactRemoved'))
                    s = projections.artifactRemoved;
                  else
                    clear s;
                    s.times = times;
                  end
                end
              else
                thisNorm = norms.artifactPresent;
                if(exist('projections', 'var'))
                  if(isfield(projections, 'artifactPresent'))
                    s = projections.artifactPresent;
                  else
                    clear s;
                  end
                end
              end
              fprintf('\nextrapolating Z scores file %d of %d, (%s)', i, length(fitFiles), char(datetime));
              jCounter = 1;
              for j = 1:length(thisNorm.means)
                fprintf('\nextrapolating Z scores file %d of %d, measure %d of %d (%s)', i, length(fitFiles), j, length(thisNorm.means), char(datetime));
                targetLabel = thisNorm.measureLabels{j};
                fitInd = find(strcmp(fit.labels, targetLabel));
                if(length(fitInd) > 0)
                  coeff = coeffvalues(fit.fits(fitInd).fun);
                  if(i == 1)
                    s.intersubjectStddev(jCounter) = thisNorm.intersubjectStddev(j);
                    s.measureLabel{jCounter} = targetLabel;
                  end
                  s.intrasubjectStddev(i, jCounter) = thisNorm.intrasubjectStddev(j);
                  %s.intrasubjectStddev(i, jCounter) = norms.rawData.artifactRemovedStddevs(i, fitInd);
                  syms x;
                  %         fprintf('\nfile %d of %d, measure %d of %d, (%s)', i, length(fitFiles), j, length(thisNorm.means), char(datetime));
                  for k = 1:kCount
                    if(solveForY)
                      if(~isnan(thisNorm.intersubjectStddev(j)))
                        y = zSamples(k) * thisNorm.intersubjectStddev(j);
                        eqn = y == coeff(1) * exp(x * coeff(2)) + coeff(3) * exp(x * coeff(4));
                        s.times(i,jCounter,k) = double(solve(eqn, x));
                      else
                        s.times(i,jCounter,k) = NaN;
                      end
                      s.hitMeasureLabel{jCounter} = targetLabel;
                    else
                      x = times(k);
                      y = coeff(1) * exp(x * coeff(2)) + coeff(3) * exp(x * coeff(4));
                      s.zValues(i,jCounter,k) = y;
                    end
                  end
                  jCounter = jCounter + 1;
                end
              end
              if(strcmp(fitFolders{folderCounter}, '/home/data/EEG/processed/Oregon/reliabilityFits/artifactRemovedReliability'))
                projections.artifactRemoved = s;
              else
                projections.artifactPresent = s;
              end
            end
          end
          save(projectedZScores, 'projections');
        else
          if(~exist('projections', 'var'))
            load(workingFile);
          end
          s = projections.artifactPresent;
        end
        
        % tabulate time needed to converge to within x standard deviations of
        % population means
        if(~solveForY)
          minTimes = NaN(size(s.zValues, 1), size(s.zValues, 2));
          for k = 1:length(targetStd)
            if(removeEyeblink)
              s = projections.artifactRemoved;
            else
              s = projections.artifactPresent;
            end
            for i = 1:size(s.zValues, 1)
              fprintf('\nfile %d of %d (%s)', i, size(s.zValues, 1), char(datetime));
              for j = 1:size(s.zValues, 2)
                z = squeeze(s.zValues(i,j,:));
                t = min(find(z < targetStd(k)));
                if(length(t) > 0)
                  minTimes(i,j,k) = s.times(t) / 60;
                else
                  minTimes(i,j,k) = max(s.times) / 60;
                end
              end
            end
            figure;
            labels = projections.artifactRemoved.measureLabel;
            isCoh = cellfun(@length, strfind(labels, 'coh')) > 0;
            isAsy = cellfun(@length, strfind(labels, 'asy')) > 0;
            isPow = cellfun(@length, strfind(labels, 'abs')) > 0;
            
            
            
            %look at asymmetry
            asyLabels = labels(isAsy);
            asyCounter = 1;
            clear asyInd asyHitLabel
            for i = 1:length(asyLabels)
              items = strsplit(asyLabels{i}, ' ');
              electrodes = strsplit(items{2}, '-');
              e1 = electrodes{1};
              e1a = e1(1:end-1);
              e1n = e1(end);
              e2 = electrodes{2};
              e2a = e2(1:end-1);
              e2n = e2(end);
              if(strcmp(e1a, e2a))
                if(e2n - e1n == 1)
                  asyInd(asyCounter) = i;
                  asyHitLabel{asyCounter} = asyLabels{i};
                  asyCounter = asyCounter + 1;
                end
              end
            end
            translation = find(isAsy);
            translation = translation(asyInd);
            asyInd = translation;
            newIsAsy = false(size(isAsy));
            newIsAsy(asyInd) = 1;
            
            
            mat = minTimes(:,isCoh | isPow,k);
            array = reshape(mat, [1, numel(mat)]);
            
            %asymmetry histogram
            if(false)
              hist(array, 50);
              if(removeEyeblink)
                eyeStatus = '(eyeblink removed)';
              else
                eyeStatus = '(eyeblink present)';
              end
              title(sprintf('Time required to converge to Z +/- %0.2f with 95%% probability %s\n', targetStd(k), eyeStatus));
              xlabel('time (minutes)');
              ylabel('number of qEEG measures');
              a = array;
              pctle = 50;
              fprintf('\nminutes for %d%% of measures to reach %f: %f\n', pctle, targetStd(k), prctile(a(~isnan(a)), pctle));
              pctle = 95;
              fprintf('\nminutes for %d%% of measures to reach %f: %f\n', pctle, targetStd(k), prctile(a(~isnan(a)), pctle));
              set(gca, 'xticklabel', {'0', '20', '40', '60', '80', '100', '>120'});
            end
            print(sprintf('/home/gmay/Documents/reliabilityFigures/fig%d', figCounter), imageType);
            figCounter = figCounter + 1;
            
            
            
            
            asyTimes = minTimes(:, asyInd, k);
            remove = any(isnan(asyTimes), 2);
            asyTimes(remove,:) = [];
            meanAsyTimes = squeeze(mean(asyTimes, 1));
            asyLabels = labels(asyInd);
            if(k==printIndex)
              freqs = {'1Hz-4Hz', '5Hz-8Hz', '9Hz-12Hz', '13Hz-24Hz', '25Hz-30Hz'};
              clear asyPlot;
              chanCounter = 1;
              for i = 1:length(freqs)
                clear chanlocs;
                chanCounter = 1;
                toPlot = find(cellfun(@length, strfind(asyLabels, freqs{i})) > 0);
                for j = 1:length(toPlot)
                  items = strsplit(asyLabels{toPlot(j)}, ' ');
                  electrodes = strsplit(items{2}, '-');
                  chanlocs(chanCounter).labels = electrodes{1};
                  chanlocs(chanCounter + 1).labels = electrodes{2};
                  asyPlot(chanCounter) = meanAsyTimes(toPlot(j));
                  asyPlot(chanCounter+1) = meanAsyTimes(toPlot(j));
                  chanCounter = chanCounter + 2;
                end
                eeg.chanlocs = chanlocs;
                eeg = setStandardLocations(eeg);
                eegChanlocs(i,:) = eeg.chanlocs;
                allAsyPlot(i,:) = asyPlot;
              end
              maxValue = max(max(allAsyPlot));
              minValue = min(min(allAsyPlot));
              fprintf('\n asymmetry: min %f max %f', minValue, maxValue);
              minValue = min(min(allAsyPlot));
              for i = 1:size(allAsyPlot, 1)
                figure;
                topoplot(allAsyPlot(i,:), eegChanlocs(i,:), 'maplimits', [minValue, maxValue]);
                colorbar;
                title(sprintf('asymmetry %s', freqs{i}));
                print(sprintf('/home/gmay/Documents/reliabilityFigures/fig%d', figCounter), imageType);
                figCounter = figCounter + 1;
              end
              
              
            end
            
          end
        end
        
        
        %plot which ones converge quickly and which ones don't.
        meanTimes = nanmean(minTimes(:, :, printIndex), 1);
        figs = plotCoherencePca(meanTimes, labels);
        for i = 1:length(figs)
          print(figs{i}, sprintf('/home/gmay/Documents/reliabilityFigures/fig%d', figCounter), imageType);
          figCounter = figCounter+1;
        end
        
        if(removeEyeblink)
          fprintf('\neyeblink artifact removed\n');
        else
          fprintf('\neyeblink artifact present\n');
        end
        
        removeEyeblink
        printIndex
        figCounter
      end
    end
  end
end
    
    % %try to simplify estimation based on ratio between intra- and
    % %inter-subject standard deviation
    % data = projections.artifactRemoved;
    % subjectCounter = 1;
    % measureCounter = 1;
    % timeCounter = 1;
    % stdRatio = data.intrasubjectStddev(subjectCounter,:) ./ data.intersubjectStddev;
    % zs = data.zValues(subjectCounter, :, timeCounter);
    % close all;
    % figure;
    % hold on;
    % plot(stdRatio);
    % plot(zs);
    % legend({'stdDev ratio', 'z scores'});
    
    
    
    
    
    %fit = load(fullfile(fitFolders{folderCounter}, fitFiles(i).name));
    
    %
    % syms x;
    %
    %
    %
    % timeEstimate.stdFraction = 0.01;
    % timeEstimate.durations = NaN(length(inputFiles), length(data.fits));
    % intermediateFolder = '/Users/Geoff/Documents/reliability5/surface digest 2/old';
    % for i = 1:length(inputFiles)
    %   data = load(fullfile(intermediateFolder, inputFiles(i).name));
    %   for j = 1:length(data.fits)
    %     fprintf('\nfile %d of %d, measure %d of %d', i, length(inputFiles), j, length(data.fits));
    %     coeffs = coeffvalues(data.fits(j).fun);
    %     fasterFirst = abs(coeffs(2)) > abs(coeffs(4));
    %     if(~fasterFirst)
    %       temp = coeffs(1:2);
    %       coeffs(1:2) = coeffs(3:4);
    %       coeffs(3:4) = temp;
    %     end
    %     coeffMat(i, j, :) = coeffs;
    %     a = coeffs(1);
    %     b = coeffs(2);
    %     c = coeffs(3);
    %     d = coeffs(4);
    %     y = norms.summary.intersubjectStd(j) * timeEstimate.stdFraction;
    %     eqn = a*exp(x*b) + c*exp(x*d) == y;
    %     sol = solve(eqn);
    %     timeEstimate.durations(i,j) = double(sol);
    %   end
    % end
    % save('exponentialFit.mat', 'coeffMat');
    % % equation is a*e^(b*x) + c*e^(d*x)
    % %   'coh icaComp1-icaComp19 9Hz-12Hz' is frickin enormous for a and c
    % %   (9e12 and -9e12)
    % maxVal = max(max(max(coeffMat)));
    % minVal = min(min(min(coeffMat)));
    % maxInd = find(coeffMat == maxVal);
    % minInd = find(coeffMat == minVal);
    % minVal = min(min(min(coeffMat)));
    %
    % extreme = abs(coeffMat) > .5;
    % coeffMat(extreme)  = 0;
    % cohMatA = squeeze(coeffMat(:,isCoh,1));
    % powMatA = squeeze(coeffMat(:,~isCoh,1));
    % cohMatB = squeeze(coeffMat(:,isCoh,2));
    % powMatB = squeeze(coeffMat(:,~isCoh,2));
    % cohMatC = squeeze(coeffMat(:,isCoh,3));
    % powMatC = squeeze(coeffMat(:,~isCoh,3));
    % cohMatD = squeeze(coeffMat(:,isCoh,4));
    % powMatD = squeeze(coeffMat(:,~isCoh,4));
    %
    % imagesc(cohMatB);
    % colorbar;
    %
    % meanD = mean(cohMatD, 2);
    % [~, indD] = sort(meanD);
    % sortCohMatD = cohMatD(indD,:);
    % imagesc(sortCohMatD);
    %
    % if(false)
    %   %find extremes
    %   %here is a very low d (low non-stationarity?)
    %   minD = cohMatD == min(min(cohMatD));
    %   minD1 = find(any(minD,2));
    %   minD2 = find(any(minD,1));
    %   i = minD1;
    %   ind = minD2;
    %   a = load(fullfile(folder, files(i).name));
    %   b = a.plot95Raw{ind};
    %
    %   %here is a very low b (low measurement noise?)
    %   minB = cohMatB == min(min(cohMatB));
    %
    %   minB1 = find(any(minB,2));
    %   minB2 = find(any(minB,1));
    %   i = minB1;
    %   ind = minB2;
    %   a = load(fullfile(folder, files(i).name));
    %   b = a.plot95Raw{ind};
    %
    %
    % end
    %
    % %check for neuropsych correlations
    % neuropsychDataFilename = '~/Documents/MATLAB/EEG/PTSD MIND for Geoff.mat';
    % neuropsychData = load(neuropsychDataFilename);
    % neuropsychData = neuropsychData.neuropsychData;
    %
    % for i = 1:length(inputFiles)
    %   [inputFolder, file, ext] = fileparts(inputFiles(i).name);
    %   file = str2num(file(4:5));
    %   rows(i) = find(neuropsychData{:,1} == file);
    % end
    % for i = 1:size(neuropsychData, 2)
    %   x = meanD;
    %   y = neuropsychData{rows, i};
    %   [rhos(i), ps(i)] = corr(x,y);
    % end
    % pInd = ps < 0.05
    %
    %
    % hist(cohMat);
    % legend({'a', 'b', 'c', 'd'});
    
    
    
    
