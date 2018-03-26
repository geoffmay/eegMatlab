clear;

plotResult = 0;
onlyRest = 0;
onlyEyesOpenCompleters = 0;
only003 = 0;
doDownsampling = 1;
concatenated = 1;
saveResult = 1;
reduceDiskSize = 0;
reducedSampleSize = 0;
calculateVarianceExplained = 1;
doDummyFile = 0;

downSampleRate = 16;

addpath('/home/data/EEG/scripts/eeglab13_4_4b/functions/sigprocfunc/');
addpath('/home/data/EEG/scripts/eeglab13_4_4b/functions/resources');
addpath(genpath('/home/gmay/Downloads/FastICA_25'));


inputFolder = '/media/eegDrive/';
%inputFolder = '/media/eegDrive/diffRef/';
%outputFolder = '/home/data/EEG/processed/Robi/coherenceIca/rest/';
outputFolder = '/media/eegDrive/ica/';
if(~exist(outputFolder, 'dir'))
    mkdir(outputFolder);
end

inputFiles = dir(inputFolder);
inputFiles([inputFiles.isdir]) = [];

if(onlyRest)
    %filter
    outcome = cellfun(@length, strfind({inputFiles.name},'outcome'));
    baseline = cellfun(@length, strfind({inputFiles.name},'baseline'));
    remove = find(~baseline & ~outcome);
    inputFiles(remove) = [];
end
if(onlyEyesOpenCompleters)
    eyesOpen = cellfun(@length, strfind({inputFiles.name},'open'));
    remove = find(~eyesOpen);
    inputFiles(remove) = [];
end
if(only003)
    is003 = cellfun(@length, strfind({inputFiles.name},'_003'));
    remove = find(~is003);
    inputFiles(remove) = [];
end

for fileNumber = 1:length(inputFiles)
    inputFile = inputFiles(fileNumber).name;
    inputPath = fullfile(inputFolder, inputFile);
    fprintf('\n%s file# %d of %d (%s)', char(datetime), fileNumber, length(inputFiles), inputFile);
    truncated = inputFile(1:strfind(inputFile,'_630')-1);
    if(length(truncated) > 0)
        outputFile = sprintf('%s%s.mat', outputFolder, truncated);
    else
        outputFile = sprintf('%s%s', outputFolder, inputFile);
    end
    if(~concatenated)
        outputFile = strrep(outputFile, '.mat', 'Fast.mat');
    else
        if(only003)
            outputFile = fullfile(outputFolder, 'master003Ica.mat');
        else
            outputFile = fullfile(outputFolder, 'masterAllSubjects.mat');
        end
    end
    if(concatenated)
        if(exist(outputFile, 'file'))
            delete(outputFile);
        end
    end
    %     if(~exist(outputFile, 'file'))
    if(doDummyFile)
        %log activity
        dummy = sprintf('started on: %s', char(datetime));
        save(outputFile, 'dummy');
    end
    output = loadCoherenceTimecourse(inputPath);
    if(isfield(output,'allDerivations'))
        allDerivations = output.allDerivations;
        fileInfo.fileDuration = size(allDerivations, 1);
        if(~exist('cohIca','var'))
            cohIca.allLabels = output.allLabels;
        else
            for i = 1:length(output.allLabels)
                if(~strcmp(cohIca.allLabels{i}, output.allLabels{i}))
                    error('channelLabel mismatch');
                end
            end
        end
    else
        fileInfo.fileDuration = 0;
    end
    if(length(output) > 0)
        fileInfo.noise = output.noiseRatio;
    end    
    fileInfo.filename = inputPath;
    metadata(fileNumber) = fileInfo;
    %             coherenceData = load(inputPath);
    %             [baseFolder,baseFile,baseExt]=fileparts(coherenceData.filename);
    %             [hasNoise, noiseRatio] = checkForMainsNoise({[baseFile,baseExt]}, baseFolder);
    %             badRefCutoff = checkForBadReference(coherenceData.channelPairs(29).coherence);
    %             if(~hasNoise)
    %             if(isfield(coherenceData, 'timeCourse'))
    
    %             %determine output filename
    %             [folder, file, ext] = fileparts(coherenceData.timeCourse.basefile);
    %             file = strrep(folder, '/home/data/EEG/data/ROBI/', '');
    %             file = strrep(file, '/', ' ');
    %             freqLabels = {'delta','theta','alpha','beta','hibeta'};
    %             singleLabels = antChannelLocs;
    %             singleLabels = singleLabels(1:33);
    %
    %             %remove mastoids
    %             m1 = cellfun(@length, strfind(singleLabels, 'M1'));
    %             m2 = cellfun(@length, strfind(singleLabels, 'M2'));
    %             cpz = cellfun(@length, strfind(singleLabels, 'CPz'));
    %             removeSingle = find(m1 | m2 | cpz);
    %             removeDouble = zeros(1,528);
    %             doubleLabels = cell(size(removeDouble));
    %             counter = 1;
    %             for i = 1:33
    %                 badI = length(find(removeSingle==i)) > 0;
    %                 for j = i+1:33
    %                     badJ = length(find(removeSingle==j)) > 0;
    %                     if(badI | badJ)
    %                         removeDouble(counter) = 1;
    %                     end
    %                     doubleLabels{counter} = sprintf('%s-%s', singleLabels{i}, singleLabels{j});
    %                     counter = counter + 1;
    %                 end
    %             end
    %             coherenceData.timeCourse.coherencePlot(:,:,find(removeDouble)) = [];
    %             coherenceData.timeCourse.powerPlot(:,:,removeSingle) = [];
    %             doubleLabels(find(removeDouble)) = [];
    %
    %
    %             %reshape data with proper labels
    %             cohSize = size(coherenceData.timeCourse.coherencePlot);
    %             powSize = size(coherenceData.timeCourse.powerPlot);
    %             allSize = [cohSize(1), cohSize(2) * cohSize(3) + powSize(2) * powSize(3)];
    %             allDerivations = NaN(allSize);
    %             allLabels = cell(1, allSize(2));
    %             freqLabels = {'delta','theta','alpha','beta','hibeta'};
    %
    %             start = 1;
    %             finish = start + 4;
    %             for i = 1:cohSize(3)
    %                 allDerivations(:, start:finish) = coherenceData.timeCourse.coherencePlot(:,:,i);
    %                 for(j = 1:length(freqLabels))
    %                     allLabels{start + j - 1} = sprintf('%s %s', doubleLabels{i}, freqLabels{j});
    %                 end
    %                 start = start + 5;
    %                 finish = start + 4;
    %             end
    %             for i = 1:powSize(3)
    %                 allDerivations(:, start:finish) = coherenceData.timeCourse.coherencePlot(:,:,i);
    %                 for j = 1:length(freqLabels)
    %                     allLabels{start + j - 1} = sprintf('%s %s abs', singleLabels{i}, freqLabels{j});
    %                 end
    %                 start = start + 5;
    %                 finish = start + 4;
    %             end
    
    if(~concatenated)
        [cohIca.icaSig, cohIca.mixingMatrix, cohIca.separatingMatrix] = fastica(allDerivations');
        close all;
        
        %                 %get channel labels only (no frequency)
        %                 onlyLabels = allLabels;
        %                 for i = 1:length(onlyLabels)
        %                     label = onlyLabels{i};
        %                     space = strfind(label, ' ');
        %                     if(length(space) > 0)
        %                         label(space(1):end) = [];
        %                     end
        %                     onlyLabels{i} = label;
        %                 end
        %                 cohIca.labelChannelOnly = onlyLabels;
        cohIca.labelChannelFreq = allLabels;
        cohIca.filename = coherenceData.timeCourse.basefile;
        cohIca.reference = coherenceData.timeCourse.reference;
        
        %                 if(calculateVarianceExplained)
        %                     fprintf('\ncalculating explained variance...');
        %                     actual = allDerivations';
        %                     actualMeans = repmat(mean(actual, 2), [1, size(actual,2)]);
        %                     actualError = actual - actualMeans;
        %                     actualVariancePlot = sqrt(sum(actualError .* actualError, 1));
        %                     originalVariance = sum(actualVariancePlot .* actualVariancePlot) / (size(actual,2)-1);
        %                     explainedVariance = NaN(1, size(cohIca.icaSig,1));
        %                     tic;
        %                     for componentIndex = 1:size(cohIca.icaSig, 1)
        %                         fprintf('\ncalculating explained variance %d of %d (%s)', componentIndex, size(cohIca.icaSig,1), char(datetime));
        %                         %icaSig = 2175 x 84722
        %                         %mixing = 2325 x 2175
        %                         %separating = 2175 x 2325
        %                         %allDerivations' = 2325 x 84722
        %                         explained = cohIca.mixingMatrix(:, componentIndex) * cohIca.icaSig(componentIndex, :);
        %                         residual = actualError - explained;
        %                         residualMeans = repmat(mean(residual, 2), [1, size(residual,2)]);
        %                         residualError = residual - actualMeans;
        %                         residualVariancePlot = sqrt(sum(residualError .* residualError, 1));
        %                         newVariance = sum(residualVariancePlot .* residualVariancePlot) / (size(actual,2)-1);
        %                         explainedVariance(componentIndex) = (originalVariance - newVariance) / originalVariance;
        %                         if(false)
        %                             %debug
        %                             close all;
        %                             figure;
        %                             hold on;
        %                             plot(residualVariancePlot, 'b')
        %                             plot(actualVariancePlot, 'r')
        %                             drawnow;
        %                             %end debug
        %                         end
        %                         fprintf(': %.1f%%', explainedVariance(componentIndex) * 100);
        %                     end
        %                     toc;
        %                     cohIca.explainedVariance = explainedVariance;
        %                 end
        
        if(saveResult)
            save(outputFile, 'cohIca', 'allLabels', '-v7.3');
        end
    else  %concatenate
        close all;
        if(length(allDerivations) > 0)
            allDerivations1 = downsample(allDerivations,downSampleRate);
            if(~exist('masterData','var'))
                masterData = allDerivations1';
            else
                masterData(:, end+1:end+size(allDerivations1,1)) = allDerivations1';
            end
        end
    end
    
    %             end
    %             end
    %     end
end

if(concatenated)
    clear allDerivations;
    [cohIca.icaSig, cohIca.mixingMatrix, cohIca.separatingMatrix] = fastica(masterData);
    cohIca.labelChannelFreq = allLabels;
    cohIca.baseFileInfo = metadata;
    cohIca.reference = 'mastoids';
    
    save(outputFile, 'cohIca', '-v7.3');
end
