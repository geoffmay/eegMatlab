inputFolder = '/Users/Geoff/Documents/avgValues';
outputFolder = '/Users/Geoff/Documents/MATLAB/processed/Oregon/avgValues';
doPlot = true;
outputFilename = fullfile(outputFolder, 'subjectNorms.mat');


%fetch means and standard deviations from all files
intermediateFilename = fullfile(outputFolder, 'intermediate.mat');
if(~exist(intermediateFilename, 'file'))
files = dir(inputFolder);
files([files.isdir]) = [];
fileCounter = 1;
while(fileCounter < length(files))
    
    filename = files(fileCounter).name;
    data = load(fullfile(inputFolder, filename));
    
    data.summary.surfSummary
    
    intermediate.meanValues(fileCounter, :) = data.summary.surfSummary.meanValue;
    intermediate.stdValues(fileCounter, :) = data.summary.surfSummary.stdValue;
    if(fileCounter == 1)
        intermediate.labels = data.summary.surfSummary.labels;
    else
        %todo: make sure labels match
    end
    intermediate.filenames{fileCounter} = filename;
    
    fileCounter = fileCounter + 1;
end
save(intermediateFilename, 'intermediate');
else
    load(intermediateFilename);
end

%compute inter and intra subject variability
summary.intersubjectMean = mean(intermediate.meanValues);
summary.intersubjectStd = std(intermediate.meanValues);
summary.intrasubjectStd = mean(intermediate.stdValues);
summary.labels = intermediate.labels;


summary.volatility = summary.intrasubjectStd ./ summary.intersubjectStd;
[summary.prunedVolatility, summary.prunedLabels] = pruneCoherenceMeasures(volatility, summary.labels);
plotCoherencePca(volatility, summary.labels)
plotCoherencePca(prunedVolatility, prunedLabels)

if(doPlot)
    figure;
    hist(intermediate.meanValues(:,1))
    
   figure;
   %    plot([summary.intersubjectMean', summary.intersubjectStd', summary.intraSubjectStd']);
   %    legend({'inter-mean', 'inter-std', 'intra-std'});
   plot([summary.intersubjectStd', summary.intrasubjectStd']);
   legend({'inter-std', 'intra-std'});
end


save(outputFilename, 'summary');
