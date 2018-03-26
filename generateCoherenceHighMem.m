clear;

%folders
logFilename = '/home/data/EEG/processed/Robi/coherencePca/log/log.mat';
outputFolder = '/media/eegDrive';

%flags
createDummyFile = 1; %"reserves" output file if multiple instances are running
addPower = 1; %add power if it's missing
doPhaseAngle = 0; 
simOnly = 0;

%constants
sampleRate = 2048;
epochLength = 1 * sampleRate;
chan1 = 100;
startIndex = (chan1-1) * epochLength + 1;
endIndex = startIndex + epochLength - 1;

hiPassHz = 1;
loPassHz = 40;

%just testing permissions
oldFolder = pwd;
cd(outputFolder);
cd(oldFolder);

% filenames = [{'/home/data/EEG/data/ROBI/ROBI_003/baseline eyes open/630158995692243270.eegData'}, ...
%   {'/home/data/EEG/data/ROBI/ROBI_003/outcome eyes open/630230539149591228.eegData'}];

%filename = '/Users/Geoff/Documents/MATLAB/EEG/Coe Collection/Robi/ROBI_003/tx 1/630165006453007385.eegData';
%        filename = '/home/data/EEG/data/ROBI/ROBI_003/outcome eyes open/630230539149591228.eegData';

filenames = getRobiDataFiles();

for fileCounter = 1:length(filenames)
    load(logFilename);
    if(~exist('log','var'))
        log.filenames = filenames;
        log.status = repmat({'pending'}, 1, length(filenames));
    end
    filename = filenames{fileCounter};
    fileIndex = find(strcmp(log.filenames, filename));
    if(length(fileIndex)==0)
        fileIndex = length(log.filenames)+1;
        log.status{fileIndex} = 'pending';
        log.filenames{fileIndex} = filename;
    end
    if(strcmp(log.status{fileIndex}, 'pending'))
        log.status{fileCounter} = sprintf('started %s', char(datetime));
        save(logFilename, 'log');
        outputStart = strfind(filename, 'ROBI_');
        if(length(outputStart) < 1)
            error('invalid filename');
        end
        outFile = filename(outputStart(1):end);
        outFile = strrep(outFile,'/','_');
        outFile = strrep(outFile,'.eegData','coherenceStats.mat');
        outputPath = fullfile(outputFolder, outFile);
        if(~exist(outputPath, 'file'))
            if(createDummyFile)
                dummy = sprintf('started on %s', char(datetime));
                save(outputPath, 'dummy');
            end
            
            channelCount = 34;
            file = dir(filename);
            fileLength = file.bytes / 8;
            fileId = fopen(filename);
            contents = fread(fileId, fileLength, 'double');
            fclose(fileId);
            
            %truncate if necessary
            sampleCount = fileLength / channelCount;
            if(sampleCount ~= floor(sampleCount))
                sampleCount = floor(sampleCount);
                fileLength = sampleCount * channelCount;
                contents = contents(1:fileLength);
            end
            labels = antChannelLocs;
            data = reshape(contents, channelCount, fileLength / channelCount)';
            clear contents;
            
            cpzIndex = find(strcmp(labels,'CPz'));
            m1Index = find(strcmp(labels,'M1'));
            m2Index = find(strcmp(labels,'M2'));
            
            maxChannel = 33;
            for i=1:size(data,1)
                sample = data(i,:);
                linkedMastoids = (sample(m1Index) + sample(m2Index)) / 2;
                avg = mean(sample(1:33));
                newSample = sample - linkedMastoids;
                avgSample = sample - avg;
                data(i,1:33) = newSample(1:33);
            end
            
            
            
            channelPair.coherence = [];
            channelPair.label = '';
            pairCounter = 1;
            
            stat.meanCoherences = NaN(maxChannel,maxChannel,5);
            stat.stdDevCoherences = NaN(maxChannel,maxChannel,5);
            stat.skewCoherences = NaN(maxChannel,maxChannel,5);
            stat.kurtosisCoherences = NaN(maxChannel,maxChannel,5);
            if(doPhaseAngle)
                stat.meanPhaseAngle = NaN(maxChannel,maxChannel,5);
                stat.stdDevPhaseAngle = NaN(maxChannel,maxChannel,5);
                stat.skewPhaseAngle = NaN(maxChannel,maxChannel,5);
                stat.kurtosisPhaseAngle = NaN(maxChannel,maxChannel,5);
            end
            stat.channelLabels = labels(1:maxChannel);
            stat.filename = filename;
            
            if(doPhaseAngle)
                [ coherencePlot, absolutePowerPlot, phaseAnglePlot ] = coherenceHighMem(data);
            else
                [ coherencePlot, absolutePowerPlot ] = coherenceHighMem(data);
            end
            
            %left off here.
            
            %             for chan1 = 1:(maxChannel)
            %                 channelPair.coherence = NaN(1,1);
            %                 channelPair.phaseAngle = NaN(1,1);
            %                 for chan2 = chan1+1:(maxChannel)
            %                     fprintf('%d, %d\n',i,j);
            %                     stat.meanCoherences(chan1,chan2,:) = mean(channelPair.coherence);
            %                     stat.stdDevCoherences(chan1,chan2,:) = std(channelPair.coherence);
            %                     stat.skewCoherences(chan1,chan2,:) = skewness(channelPair.coherence);
            %                     stat.kurtosisCoherences(chan1,chan2,:) = kurtosis(channelPair.coherence);
            %                     if(doPhaseAngle)
            %                         stat.meanPhaseAngle(chan1,chan2,:) = mean(channelPair.phaseAngle);
            %                         stat.stdDevPhaseAngle(chan1,chan2,:) = std(channelPair.phaseAngle);
            %                         stat.skewPhaseAngle(chan1,chan2,:) = skewness(channelPair.phaseAngle);
            %                         stat.kurtosisPhaseAngle(chan1,chan2,:) = kurtosis(channelPair.phaseAngle);
            %                     end
            %                     stat.coherenceThreshold = .9;
            %                     supraThreshold = channelPair.coherence > stat.coherenceThreshold ;
            %                     for i = 1:size(supraThreshold,2)
            %                         a = supraThreshold(:,i);
            %                         b = diff(a);
            %                         rising = find(b == 1);
            %                         falling = find(b == -1);
            %                         while(isFallingLessThanRising(falling,rising))
            %                             falling(1) = [];
            %                         end
            %                         while(length(rising) > length(falling))
            %                             rising(end) = [];
            %                         end
            %                         while(length(rising) < length(falling))
            %                             falling(end) = [];
            %                         end
            %                         %chop the end off because coherence has some lag in finding quality
            %                         %data.
            %                         falling = falling - 64;
            %                         if((length(rising) > 0) && (length(falling) > 0))
            %                             epochLengths = falling - rising;
            %                         else
            %                             epochLengths = [];
            %                         end
            %                         minEpochLength = 129;
            %                         goodEpoch = find(epochLengths > minEpochLength);
            %                         rising = rising .* 16 + 1024;
            %                         falling = falling .* 16 + 1024;
            %                         freqLabels = [{'Delta'}, {'Theta'}, {'Alpha'}, {'Beta'}, {'Hibeta'}];
            %                         risingLabel = sprintf('supraThreshold%sEpochStart', freqLabels{i});
            %                         fallingLabel = sprintf('supraThreshold%sEpochEnd', freqLabels{i});
            %                         stat.(risingLabel) = rising(goodEpoch);
            %                         stat.(fallingLabel) = falling(goodEpoch);
            %                     end
            %
            %
            %                     if(any(any(isnan(channelPair.coherence))));
            %                         fprintf('i = %d, j = %d\n', chan1, chan2);
            %                         error('unexpected coherenceValue');
            %                     end
            %                     channelPair.label =  sprintf('%s-%s',labels{chan1},labels{chan2});
            %                     progress = sprintf('(%s) file %d/%d, channelPair %d/%d\n',char(datetime), fileCounter, length(filenames),pairCounter,maxPairs);
            %                     fprintf(progress);
            %                     save(outputPath, 'progress');
            %                     channelPairs(pairCounter) = channelPair;
            %                     pairCounter = pairCounter + 1;
            %                     if(false)
            %                         %debug
            %                         normalizedPhaseAngle = channelPair.phaseAngle ./ (2*pi) + .5;
            %                         figure;
            %                         hold on;
            %                         plot(channelPair.coherence(:,1), 'r');
            %                         plot(normalizedPhaseAngle, 'b');
            %                         legend('coherence', 'phase');
            %                         pan xon;
            %                         zoom xon;
            %                         %end debug
            %                     end
            %                 end
            %             end
            singleLabels = antChannelLocs;
            for channelIndex = 1:33
                channels(channelIndex).label = singleLabels{channelIndex};
            end
            save(outputPath, 'channelPairs', 'channels', 'filename', '-v7.3');
        else %outputfile exists
            if(addPower)
                fprintf('\n(%s) %d of %d', char(datetime), fileCounter, length(filenames));
                load(outputPath);
                base = loadRobiDataFile(filename);
                for channelIndex = 1:33
                    [channels(channelIndex).absolutePower, channels(channelIndex).relativePower] = smearedPower(base(:,channelIndex));
                end
                if(false)
                    %plot timecourse
                    x = (.5:(1/128):.5+(size(channels(1).absolutePower,1)-1)/128)'
                    whiteStep = 0.15;
                    blue = [0 0 .5; .5 0 0];
                    colorCounter = 1;
                    for c1 = 1:size(blue,1)
                        for c = 1:6
                            for c3 = 1:3
                                blues(colorCounter,c3) = blue(c1, c3) + whiteStep * c * (1-blue(c1,c3));
                            end
                            colorCounter = colorCounter + 1;
                        end
                    end
                    figure;
                    hold on;
                    set(gca,'colororder',blues)
                    plot(x, channels(1).absolutePower)
                    rel = channels(1).relativePower;
                    maxRel = max(max(rel));
                    minRel = min(min(rel));
                    maxAbs = max(max(channels(1).absolutePower));
                    minAbs = min(min(channels(1).absolutePower));
                    ratio = (maxAbs - minAbs) / (maxRel - minRel);
                    for i = 1:size(rel,1)
                        for j = 1:size(rel,2)
                            rel(i,j) = (rel(i,j) - minRel) * ratio + minAbs;
                        end
                    end
                    
                    plot(x, rel);
                    pan xon;
                    zoom xon;
                    
                end
                save(outputPath, 'channelPairs', 'channels', 'filename', '-v7.3');
            end
        end
        load(logFilename);
        log.status{fileCounter} = sprintf('finished %s', char(datetime));
        save(logFilename, 'log');
    end
end
