clear;
clc;

doGranger = false;

%get list of files
inputFolder = '/Users/Geoff/Documents/MATLAB/EEG/rebuilt feedback'; %mac
inputFolder = '/Volumes/My Book/zScoresComplete'; %mac ext hdd
% folder = '/home/data/EEG/processed/Robi/fullZ'; %server
inputFolder = '/media/HERMANOSBAC/data/EEG1'; %server external drive

outputFolder = fullfile(inputFolder, 'pca');
outputFolder = '/home/data/EEG/processed/Robi/fullZ/pca';

files = dir(inputFolder);
files = files(~[files.isdir]);

%get variable info from one of the notes files.
searchString = '.zScorenotes.txt';
searchHits = strfind({files.name}, searchString);
searchHits = find(cellfun(@length, searchHits));
fileId = fopen(fullfile(inputFolder,files(searchHits(1)).name));
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
    removeSuffix = '.zScorenotes.txt1.chunk';
    if(length(filename) > length(removeSuffix))
        if(filename(end-length(removeSuffix)+1:end) == removeSuffix)
            remove(i) = 1;
        end
    end
end
files(find(remove)) = [];

%process all files
array = [];
oldFilename = files(1).name;
for i = 1:length(files)
    simpleName = oldFilename(1:strfind(oldFilename, '.zScore')-1);
    outputPath = fullfile(outputFolder, sprintf('%s.nipals.mat',simpleName));
    if(length(strfind(files(i).name, '.zScore1.chunk')) > 0 && length(array) > 0)

        %shape data matrix based on column headers from info file
        columns = length(columnLabels);
        rows = length(array) / columns;
        contents = reshape(array, columns, rows)';
        array = [];
        maxJ = size(contents, 2);
        
        eegMeta.sampleRateHertz = 2048;
        eegMeta.neuroguideSampleRateHertz = 128;
        
        %pca
        tic;
        fprintf('\n\ndoing pca on file %d of %d', i, length(files));
        dataPca.maxIterations = power(10,6);
        dataPca.desiredComponents = 35;
        [dataPca.T, dataPca.P, dataPca.pcvar] = nipals(...
            contents(:,1:end-8), dataPca.desiredComponents, ...
            dataPca.maxIterations);
        dataPca.elapsed = toc;
        dataPca.description = 'nipals-based primary component analysis of all neuroguide parameters';
        
        %granger
        if(doGranger) %granger causality
            tickSize = floor(maxJ/100);
            correlationRs = NaN(1,maxJ);
            correlationPs = NaN(1,maxJ);
            grangerF = NaN(maxJ, 1);
            grangerCrit = NaN(maxJ, 1);
            grangerCausal = NaN(maxJ, 1);
            tic;
            for j = 1:maxJ
                fprintf('.');
                if(mod(j,tickSize) == 0)
                    fprintf('.');
                    if(mod(j,tickSize * 10))
                        fprintf('%d', j / tickSize);
                    end
                end
                granger.max_lag = 5;
                granger.alpha = 0.05 / length(columnLabels);
%                 x = contents(:,j);
%                 y = contents(:,maxJ-3); %feedback
                y = contents(:,j);
                x = dataPca.T(:,1); %first Pca component time series
                [correlationRs(j), correlationPs(j)] = corr(x,y);
                [grangerF(j), grangerCrit(j)] = granger_cause(x,y,granger.alpha,granger.max_lag);
                grangerCausal(j) = grangerF(j) > grangerCrit(j);
            end
            
            causalRatio = grangerF ./ grangerCrit;
            standardDev = std(contents,1)';
            meanValue = mean(contents,1)';
            normalizedCausalRatio = causalRatio ./ standardDev;
            tab = table(columnLabels, grangerF, grangerCrit, grangerCausal, causalRatio, meanValue, standardDev, normalizedCausalRatio);
            granger.table = tab;
            granger.elapsed = toc;
            granger.description = 'the granger-causal effect of each neuroguide parameter to the first pca component';
            save(outputPath, 'dataPca', 'granger', 'eegMeta');
        else
            save(outputPath, 'dataPca', 'eegMeta');
            
%         else
%             contents = contents(:, 1:end-10);
%             save(outputPath, 'dataPca', 'columnLabels');
        end
        oldFilename = files(i).name;
    end
    
    peekName = files(i).name;
    peekName = peekName(1:strfind(peekName, '.zScore')-1);
    peekPath = fullfile(outputFolder, sprintf('%s.nipals.mat',peekName));
    if(~exist(peekPath, 'file'))
        %read the file
        inputFilename = fullfile(inputFolder, files(i).name);
        fileId = fopen(inputFilename);
        contents = fread(fileId, files(i).bytes/4, 'single');
        array = [array; contents];
        fclose(fileId);
    end
    %     end
    
    %     %unzip the file if necessary
    %     if(strcmp(files(i).name(end-3:end), '.zip'))
    %         tempFilename = fullfile(outputFolder, files(i).name);
    %         [catCode, catResult] = system(sprintf('cat %s* > %s', inputFilename(1:end-2), tempFilename));
    %         [zipCode, zipResult] = system(sprintf('unzip %s -d %s*', tempFilename, outputFolder));
    %         if(zipCode ~=0)
    %             error('error unzipping file %s', tempFilename);
    %         end
    %         target = 'extracting: ';
    %         index = strfind(zipResult, target) + length(target);
    %         inputFilename = strtrim(zipResult(index:end));
    %         stuffAddedByUnzip = inputFilename(length(outputFolder)+2:end);
    %     end
    %     outputPath = fullfile(outputFolder, sprintf('%s.nipals.mat',files(i).name(1:end-4)));
    
    %execute only if not already done
    
    
%     %clean up any unzipping that was done
%     if(strcmp(files(i).name(end-3:end), '.zip'))
%         delete(tempFilename);
%         slashes = strfind(stuffAddedByUnzip, '/');
%         if(length(slashes) > 0)
%             tempRoot = fullfile(outputFolder, stuffAddedByUnzip(1:slashes(1)));
%             rmdir(tempRoot, 's');
%         else
%             delete(fullfile(outputFolder, stuffAddedByUnzip));
%         end
%     end
%     elapsed = toc;
    
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


