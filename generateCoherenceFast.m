
addPower = 1;
doRereference = 0;
doPhaseAngle = 1;
simOnly = 0;
summaryOnly = 0;
doAllChannel = 1;
saveResult = 1;
onlyEyesOpenCompleters = 1;


sampleRate = 2048;
epochLength = 1 * sampleRate;
chan1 = 100;
startIndex = (chan1-1) * epochLength + 1;
endIndex = startIndex + epochLength - 1;

hiPassHz = 1;
loPassHz = 40;
filenames = getRobiDataFiles();

if(onlyEyesOpenCompleters)
    eyesOpen = cellfun(@length, strfind(filenames,'open'));
    remove = find(~eyesOpen);
    filenames(remove) = [];
end

for fileCounter = 1:length(filenames)
filename = filenames{fileCounter};
if(summaryOnly)
    outputFolder = '/home/data/EEG/processed/Robi/coherenceReref';
else
    outputFolder = '/home/data/EEG/processed/Robi/coherence';
    oldFolder = pwd;
    outputFolder = '/media/eegDrive/csd';
    outputFolder = '/media/eegDrive';
    outputFolder = '/media/eegDrive/fast';
    cd(outputFolder);  %just testing permissions
    cd(oldFolder);
end
outputStart = strfind(filename, 'ROBI_');
if(length(outputStart) < 1)
    error('invalid filename');
end
outFile = filename(outputStart(1):end);
outFile = strrep(outFile,'/','_');
chaffIndex = strfind(outFile,'_63');
outFile = outFile(1:chaffIndex-1);
outFile = sprintf('%s.mat', outFile);
outFile = strrep(outFile,'.eegData','coherenceStats.mat');

channelCount = 34;
file = dir(filename);
fileLength = file.bytes / 8;
sampleCount = fileLength / channelCount;
if(sampleCount ~= floor(sampleCount))
    sampleCount = floor(sampleCount);
end
if(~saveResult)
    %load only part of the file for speed
    sampleCount = 2048 * 30;
end
fileLength = sampleCount * channelCount;
fileId = fopen(filename);
contents = fread(fileId, fileLength, 'double');
fclose(fileId);

labels = antChannelLocs;
data = reshape(contents, channelCount, fileLength / channelCount)';
clear contents;
cpzIndex = find(strcmp(labels,'CPz'));
m1Index = find(strcmp(labels,'M1'));
m2Index = find(strcmp(labels,'M2'));


outputFolder = '/media/eegDrive/diffRef';
references = {'Cpz', 'Mastoids', 'Average'}; %todo: laplace
for referenceCounter = 1:length(references)
    reference = references{referenceCounter};
    
    outFile1 = strrep(outFile,'.mat',sprintf('-%s.mat', reference));
    outputPath = fullfile(outputFolder, outFile1);
    %if(~exist(outputPath, 'file'))
    if(true)
        
        EEG.data = data;
        EEG.data(:,34) = [];
        
        for i=1:size(EEG.data,1)
            sample = EEG.data(i,:);
            if(strcmp(reference, 'Cpz'))
                newReference = 0;
            elseif(strcmp(reference, 'Mastoids'))
                newReference = (sample(m1Index) + sample(m2Index)) / 2;
            elseif(strcmp(reference, 'Average'))
                newReference = mean(sample);
            end
            newSample = sample - newReference;
            EEG.data(i,1:33) = newSample(1:33);
        end
        
        EEG.data = EEG.data';
        EEG.srate = 2048;
        EEG.nbchan = size(EEG.data,1);
        [ timeCourse.coherencePlot, timeCourse.x, timeCourse.powerPlot ] = allChannelCoherence( EEG );
        timeCourse.reference = reference;
        timeCourse.basefile = filename;
        if(saveResult)
            save(outputPath, 'timeCourse', '-v7.3');
        else
            combinedData.derivation = timeCourse;
            combinedData.raw = EEG.data;
            combinedData.reference = reference;
            if(~exist('allData', 'var'))
                allData = combinedData;
            else
                allData(end+1) = combinedData;
            end
        end
    end
end


if(~saveResult)
    labels = antChannelLocs;
    m1 = find(strcmp(labels, 'M1'));
    m2 = find(strcmp(labels, 'M2'));
    cpz = find(strcmp(labels, 'CPz'));
    removeSingle = [m1 m2 cpz];
    removeDouble = [];
    counter = 1;
    for i = 1:33
        badI = any(i==removeSingle);
        for j = i+1:33
            badJ = any(j==removeSingle);
            if(badI | badJ)
                removeDouble(end+1) = counter;
            end
            counter = counter + 1;
        end
    end
    for i = 1:length(allData)
        allData(i).raw(removeSingle,:) = [];
        allData(i).derivation.coherencePlot(:,:,removeDouble) = [];
        allData(i).derivation.powerPlot(:,:,removeSingle) = [];
        
    end
end
end