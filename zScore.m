function [ zScores ] = zScore( dataTimeCourse, age, eyesOpen, frequencyOfChoice )
%ZSCORE Summary of this function goes here
%   Detailed explanation goes here



if(~exist('frequencyOfChoice','var'))
    frequencyOfChoice = 3;
end
if(~exist('age','var'))
    age = 30;
end
if(~exist('eyesOpen','var'))
    eyesOpen = 1;
end
if(eyesOpen)
    eyes = 'open';
else
    eyes = 'closed';
end

%load the data
filename= '/home/data/EEG/processed/Robi/zScore.mat';
load(filename);
frequencyLookupLabels = [{'Delta 1-4Hz'}, {'Theta 4-8Hz'}, {'Alpha 8-12Hz'}, {'Beta 12-25Hz'}, {'HiBeta 25-30Hz'},{'Alpha1 8-10Hz'},{'Alpha2 10-12Hz'},{'Beta1 12-15Hz'},{'Beta2 15-18Hz'},{'Beta3 18-25Hz'}];
frequencyLookupLabels = frequencyLookupLabels(1:5);
aniChanLocs = {'Fp1','Fp2','F3','F4','C3','C4','P3','P4','O1','O2','F7','F8','T3','T4','T5','T6','Fz','Cz','Pz'};
antChnLocs = antChannelLocs;

%determine rows of interest
eyeOpen = strcmp(data{:,'Subject_Eyes'}, eyes);
absPower = strcmp(data{:,'Analytic_Type'}, 'absolute power mean');
absPowerStd = strcmp(data{:,'Analytic_Type'}, 'absolute power standard deviation');
relPower = strcmp(data{:,'Analytic_Type'}, 'relative power mean');
relPowerStd = strcmp(data{:,'Analytic_Type'}, 'relative power standard deviation');
coherence = strcmp(data{:,'Analytic_Type'}, 'coherence mean');
coherenceStd = strcmp(data{:,'Analytic_Type'}, 'coherence standard deviation');

zScores.channelPairLabels = {dataTimeCourse.coherencePlot.label};
% fprintf('\nallocating... ');
cohDim = [size(dataTimeCourse.coherencePlot(1).coherence,1), 5, 171];
powDim = [size(dataTimeCourse.powerPlot(1).absolutePower,1), 5, 19];
zScores.absPower = NaN(powDim);
zScores.relPower = NaN(powDim);
zScores.coherence = NaN(cohDim);
ageIndex = 7;
foundAge = false;
while(~foundAge)
    columnLabel = data.Properties.VariableNames{ageIndex};
    minAge = str2num(columnLabel(5:end));
    if(age < minAge)
        foundAge = true;
    else
        ageIndex = ageIndex + 1;
    end
end

%%parse tables
freqRows = strcmp(data{:,'Spectral_Bin_Label'}, frequencyLookupLabels{frequencyOfChoice});
meanAbsPowRows = find(eyeOpen & absPower & freqRows);
stdAbsPowRows = find(eyeOpen & absPowerStd & freqRows);
meanRelPowRows = find(eyeOpen & relPower & freqRows);
stdRelPowRows = find(eyeOpen & relPowerStd & freqRows);
meanCohRows = find(eyeOpen & coherence & freqRows);
stdCohRows = find(eyeOpen & coherenceStd & freqRows);

absPowMeans = data{meanAbsPowRows, ageIndex};
absPowStds = data{stdAbsPowRows, ageIndex};
relPowMeans = data{meanRelPowRows, ageIndex};
relPowStds = data{stdRelPowRows, ageIndex};
cohMeans = data{meanCohRows, ageIndex};
cohStds = data{stdCohRows, ageIndex};

%%prepare labels
channelPairLabels = data{meanCohRows, 'Channel_Pair_Label'};
channelLabels = data{meanAbsPowRows, 'Channel_Pair_Label'};
if(frequencyOfChoice == 1)
    zScores.channelLabels = channelLabels;
    zScores.channelPairLabels = channelPairLabels;
end
searchLabels = channelLabels;
for i = 1:length(searchLabels)
    searchLabels = strrep(searchLabels, 'T3', 'T7');
    searchLabels = strrep(searchLabels, 'T4', 'T8');
    searchLabels = strrep(searchLabels, 'T5', 'P7');
    searchLabels = strrep(searchLabels, 'T6', 'P8');
end
channelIndexes = NaN(1, length(searchLabels));
channelPairIndexes = NaN(1, length(channelPairLabels));
for i = 1:length(channelIndexes)
    channelIndexes(i) = find(strcmp(antChnLocs, searchLabels{i}));
end
pairCounter = 1;
for i = 1:length(searchLabels)
    for j = i+1:length(searchLabels)
        pairLabel = sprintf('%s-%s', searchLabels{i}, searchLabels{j});
        index = find(strcmp({dataTimeCourse.coherencePlot.label}, pairLabel));
        if(length(index) == 0)
            pairLabel = sprintf('%s-%s', searchLabels{j}, searchLabels{i});
            index = find(strcmp({dataTimeCourse.coherencePlot.label}, pairLabel));
        end
        channelPairIndexes(pairCounter) = index;
        pairCounter = pairCounter + 1;
    end
end



totalPower = NaN(size(dataTimeCourse.powerPlot(1).absolutePower,1), length(dataTimeCourse.powerPlot));
for chanCounter = 1:length(channelIndexes)
    for timeCounter = 1:size(dataTimeCourse.powerPlot(chanCounter).absolutePower,1)
        %         fprintf('\n%d %d', chanCounter, timeCounter);
        totalPower(timeCounter, chanCounter) = sum(dataTimeCourse.powerPlot(channelIndexes(chanCounter)).absolutePower(timeCounter, :));
    end
end

for frequencyOfChoice = 1:length(frequencyLookupLabels)
    
    %%compute Z scores
    timeLength = size(dataTimeCourse.coherencePlot(1).coherence, 1);
    cohData = NaN(length(channelPairIndexes), 1);
    absPowData = NaN(length(channelIndexes), 1);
    relPowData = NaN(length(channelIndexes), 1);
    
    for timeCounter = 1:timeLength
        for chanPairCounter = 1:length(channelPairIndexes)
            cohData(chanPairCounter) = dataTimeCourse.coherencePlot(channelPairIndexes(chanPairCounter)).coherence(timeCounter, frequencyOfChoice);
        end
        
        for chanCounter = 1:length(channelIndexes)
            absPowData(chanCounter) = dataTimeCourse.powerPlot(channelIndexes(chanCounter)).absolutePower(timeCounter, frequencyOfChoice);
            relPowData(chanCounter) = absPowData(chanCounter) / totalPower(timeCounter, chanCounter);
        end
        cohData = cohData .* 100;
        cohZ = (cohData - cohMeans) ./ cohStds;
        absPowZ = (absPowData - absPowMeans) ./ absPowStds;
        relPowZ = (relPowData - relPowMeans) ./ relPowStds;
        zScores.absPower(timeCounter, frequencyOfChoice, :) = absPowZ;
        zScores.relPower(timeCounter, frequencyOfChoice, :) = relPowZ;
        zScores.coherence(timeCounter, frequencyOfChoice, :) = cohZ;
    end
end


end

