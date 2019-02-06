%filename = 'C:\Users\Neuro\Documents\MATLAB\processed\Oregon\artifactRemovedReliability\PM101Surface.mat';
function output = loadReliabilityMatrix(filename, measureLabels)

fileData = load(filename);
if(~exist('measureLabels', 'var'))
    measureLabels = fileData.summary.surfLabels;
end
%numberOfMeasures = length(fileData.summary.surfResample(1).averageDifference);
numberOfMeasures = length(measureLabels);
output.measureLabels = measureLabels;
output.averageDifference = NaN(numberOfMeasures, length(fileData.summary.surfResample));
output.percentile95Difference = NaN(numberOfMeasures, length(fileData.summary.surfResample));
stdDev95 = norminv(0.95);

for measureIndex = 1:numberOfMeasures
    %fprintf('\nmeasure %d of %d', removedMeasureIndex, numberOfMeasures);
    index = find(strcmp(fileData.summary.surfLabels, measureLabels{measureIndex}));
    %measureLabel = fileData.summary.surfLabels{measureIndex};
    %output.measureLabels{measureIndex} = measureLabel;
    for durationIndex = 1:length(fileData.summary.surfResample)
        resampled = fileData.summary.surfResample(durationIndex);
        avgDiff = resampled.averageDifference(index);
        stdDevDiff = resampled.stddevDifference(index);
        output.averageDifference(measureIndex, durationIndex) = avgDiff;
        output.percentile95Difference(measureIndex, durationIndex) = stdDev95 * stdDevDiff + avgDiff;
        if(measureIndex == 1)
            output.sampleMinuteDurations(durationIndex) = resampled.frameCount / 128 / 60;
        end
    end
end